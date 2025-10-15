
using LegendHDF5IO
using StatsBase
using ProgressMeter
using StatsBase
using CalibrationTemplateFits
using Glob

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
function read_mc_files(det::Union{Symbol,String}, files::Vector)

    outputs = nothing
    if length(files)==0
        throw(ArgumentError("No files"))
    end

    for fi in files
        out_tmp = lh5open(fi, "r") do f
            f["hit/$det/energy"][:]
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
    pattern::Union{Regex,Nothing},
)

    if (pattern == nothing)
        pattern = r".*z-offset_([-\d.]+)_phi-offset_([-\d.]+)"
    end

    n_prim = nothing
    outputs = Dict()

    @showprogress for det in dets
        histograms = HistogramWithPars[]

        zs = Set([])
        φs = Set([])

        for folder in folders

            files = glob("*", folder)
            p_name = split(folder, "/")[end]
            energies = read_mc_files(det, files)

            z, φ = extract_mc_coords(String(p_name), pattern)

            h = append!(Histogram(binning), energies)
            push!(histograms, HistogramWithPars(h, z = z, φ = φ))

            push!(zs, z)
            push!(φs, φ)

        end
        outputs[det] =
            GeneralisedHistogram(histograms, z = _set2range(zs), φ = _set2range(φs))
    end
    return outputs
end
