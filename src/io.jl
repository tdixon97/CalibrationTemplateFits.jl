
using LegendHDF5IO
using StatsBase
using ProgressMeter
using StatsBase
using CalibrationTemplateFits
using Glob
using Base.Threads
using LegendHDF5IO
"""
    function read_data_histograms(
        data_path::String,
        file_list::AbstractVector{String},
        det_map::Dict,
        binning::Union{AbstractRange,AbstractVector},
    )

Read the data from `data_path` and the files in `file_list` and make histograms with the `binning` for each detector in the keys of `dets` (with the values
being the "rawids").
"""
function read_data_histograms(
    data_path::String,
    file_list::AbstractVector{String},
    det_map::Dict,
    binning::Union{AbstractRange,AbstractVector},
)

    data = read_data(data_path, file_list)

    data_sel = data[(.!data.coincident.puls) .& (.!data.trigger.is_forced)]

    data_hists = Dict(
        det => get_data_histogram(det_map[det], data_sel, binning) for det in keys(det_map)
    )
    return data_hists

end
""""
    read_data(path::String, files::Vector)

Read the data from each file in the `files` list, at the path `path`.
"""
function read_data(path::String, files::Vector)

    data = nothing
    @showprogress for file in files

        table_tmp = lh5open("$path/$file", "r") do f
            f["evt"][:]
        end

        # add it to the rest
        if (data != nothing)
            data = vcat(data, table_tmp)
        else
            data = table_tmp
        end

    end
    return data
end
"""
    read_mc_files(det::Union{Symbol,String}, files::Vector)

Read the hit tier data from each file in `files` and the detector `det`.
"""
function read_mc_files(det::Union{Symbol,String}, files::Vector; field::String = "energy")

    outputs = nothing
    if length(files)==0
        throw(ArgumentError("No files"))
    end

    for fi in files
        out_tmp = lh5open(fi, "r") do f
            f["hit/$det/$field"][:]
        end
        if (outputs!=nothing)
            outputs = vcat(outputs, out_tmp)
        else
            outputs = out_tmp
        end
    end
    outputs
end


"""
    extract_mc_coords(folder::String,pattern::String

Extract the information from the mc file path
"""
function extract_mc_coords(folder::String, pattern::Regex)
    m = match(pattern, folder)
    if m !== nothing
        x = parse(Float64, m.captures[1])
        y = parse(Float64, m.captures[2])

    else
        throw(ArgumentError("No match found in $folder to $pattern"))
    end
    return x, y
end

"""
    read_models(dets::AbstractVector, mc_path::String, label::String,binning::Union{AbstractVector,AbstractRange})

Read the MC files from `mc_path`, saving the results to the `GeneralisedHistogram`` models. We expect every
MC file to be stored in a folder defining a 2D grid. The coordinates are extracted from the file name based on
the supplied `pattern` regex.
"""
function read_models(
    dets::AbstractVector,
    folders::AbstractVector{String},
    binning::Union{AbstractVector,AbstractRange},
    pattern::Union{Regex,Nothing} = nothing;
    vary_fccd::Bool = false,
)

    if pattern === nothing
        pattern = r".*z-offset_([-\d.]+)_phi-offset_([-\d.]+)"
    end

    n_dets = length(dets)
    outputs_array = Vector{Any}(undef, n_dets)

    @showprogress @threads for i = 1:n_dets
        det = dets[i]

        histograms = HistogramWithPars[]
        zs = Set{Float64}()
        φs = Set{Float64}()

        isempty(folders) &&
            throw(ValueError("no folders found in the input path - check this!"))

        for folder in folders
            files = glob("*", folder)
            p_name = split(folder, "/")[end]
            z, φ = extract_mc_coords(String(p_name), pattern)

            if !vary_fccd
                energies = read_mc_files(det, files; field = "energy")
                h = append!(Histogram(binning), energies)
                push!(histograms, HistogramWithPars(h, z = z, φ = φ))
                push!(zs, z)
                push!(φs, φ)
            else
                energies = read_mc_files(det, files; field = "surface/energies")
                dist = read_mc_files(det, files; field = "surface/dist_to_nplus")
                bulk_energy = read_mc_files(det, files; field = "bulk_energy")

                for fccd = 0:0.2:2
                    energy = copy(bulk_energy)
                    for idx in eachindex(energy)
                        act = CalibrationTemplateFits.piecewise_linear_activeness.(
                            Float64.(dist[idx]);
                            fccd = fccd,
                            dlf = 0.5,
                        )
                        if !isempty(act)
                            energy[idx] += sum(act .* energies[idx])
                        end
                    end
                    h = append!(Histogram(binning), energy)
                    d = Dict(Symbol("fccd_$det") => fccd)
                    kw = (z = z, φ = φ, (; d...)...)
                    push!(histograms, HistogramWithPars(h, kw))
                    push!(zs, z)
                    push!(φs, φ)
                end
            end
        end

        # build per-detector histogram
        outputs_array[i] = if !vary_fccd
            GeneralisedHistogram(histograms; z = _set2range(zs), φ = _set2range(φs))
        else
            d = Dict(Symbol("fccd_$det") => 0:0.2:2)
            kw = (z = _set2range(zs), φ = _set2range(φs), (; d...)...)
            GeneralisedHistogram(histograms, kw)
        end
    end

    # convert back to Dict keyed by det
    outputs = Dict(dets[i] => outputs_array[i] for i = 1:n_dets)

    outputs
end

""""
    read_models_evt(dets::AbstractVector, folders::AbstractVector{String}, binning::Union{AbstractVector,AbstractRange}, pattern::Union{Regex,Nothing} = nothing)

Read the models from the evt tier simulations
"""
function read_models_evt(
    dets::AbstractVector,
    folders::AbstractVector{String},
    binning::Union{AbstractVector,AbstractRange},
    det_map::Dict,
    pattern::Union{Regex,Nothing} = nothing,
)

    if pattern === nothing
        pattern = r".*z-offset_([-\d.]+)_phi-offset_([-\d.]+)"
    end

    n_dets = length(dets)
    outputs_array = Vector{Any}(undef, n_dets)
    det_histograms = Dict(det=>HistogramWithPars[] for det in dets)
    det_zs = Dict(det=>Set{Float64}() for det in dets)
    det_φs = Dict(det=>Set{Float64}() for det in dets)

    isempty(folders) &&
        throw(ValueError("no folders found in the input path - check this!"))

    #first loop over folders

    for folder in folders
        @info "   ... reading $folder"
        # read the data
        file_list = glob("*", folder)
        data = read_data("/", file_list)

        p_name = split(folder, "/")[end]
        z, φ = extract_mc_coords(String(p_name), pattern)

        for i = 1:n_dets
            det = dets[i]

            r = flatview(data.channel)
            energy = flatview(data.energy)
            is_good = flatview(data.is_good)
            energy_sel = energy[(r .== det_map[det]) .& (is_good)]

            h = append!(Histogram(binning), energy_sel)

            push!(det_histograms[det], HistogramWithPars(h, z = z, φ = φ))
            push!(det_zs[det], z)
            push!(det_φs[det], φ)

        end

    end

    outputs = Dict()
    for det in dets
        # build per-detector histogram
        outputs[det] = GeneralisedHistogram(
            det_histograms[det];
            z = _set2range(det_zs[det]),
            φ = _set2range(det_φs[det]),
        )
    end

    outputs
end

"""
     read_hist(group::String,file::String)
"""
function read_hist(group::String, file::String)

    binning = lh5open("$file", "r") do f
        f["$group/binning/axis_0"]
    end
    weights = lh5open("$file", "r") do f
        f["$group/weights"][:]
    end
    h = LegendHDF5IO._nt_to_histogram((
        binning = [binning],
        weights = weights,
        isdensity = false,
    ))

    return h

end

"""
    read_models_hist(
        dets::AbstractVector,
        files::AbstractVector{String},
        binning::Union{AbstractVector,AbstractRange},
        pattern::Union{Regex,Nothing} = nothing,
        group::String
    )

Read precomputed histograms from `files` and rebin them to the desired `binning``.
"""
function read_models_hist(
    dets::AbstractVector,
    files::AbstractVector{String},
    binning::Union{AbstractVector,AbstractRange},
    pattern::Union{Regex,Nothing} = nothing,
    group::String = "hist",
)

    if pattern === nothing
        pattern = r".*z-offset_([-\d.]+)_phi-offset_([-\d.]+)"
    end

    n_dets = length(dets)
    det_histograms = Dict(det => HistogramWithPars[] for det in dets)
    det_zs = Dict(det => Set{Float64}() for det in dets)
    det_φs = Dict(det => Set{Float64}() for det in dets)

    isempty(files) && throw(ArgumentError("no files found in the input path - check this!"))

    for file in files

        # list histogram files (one per detector, typically)

        p_name = split(split(file, "/")[end], ".")[1]
        z, φ = extract_mc_coords(String(p_name), pattern)

        for det in dets

            h = read_hist("$group/$det", file)
            h = CalibrationTemplateFits.rebin_integer(h, binning)
            push!(det_histograms[det], HistogramWithPars(h, z = z, φ = φ))

            push!(det_zs[det], z)
            push!(det_φs[det], φ)
        end
    end

    outputs = Dict()
    for det in dets
        outputs[det] = GeneralisedHistogram(
            det_histograms[det];
            z = _set2range(det_zs[det]),
            φ = _set2range(det_φs[det]),
        )
    end

    outputs
end
