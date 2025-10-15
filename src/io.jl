
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
    read_models(dets::AbstractVector, mc_path::String, label::String,binning::Union{AbstractVector,AbstractRange})

Read the MC files from `mc_path`, saving the results to the `GeneralisedHistogram`` models.
"""
function read_models(
    dets::AbstractVector,
    folders::AbstractVector{String},
    binning::Union{AbstractVector,AbstractRange},
)

    n_prim = nothing
    outputs = Dict()

    @showprogress for det in dets
        histograms = HistogramWithPars[]

        zs = Set([])
        φs = Set([])

        for folder in folders

            files = glob("*", folder)
            p_name = split(folder, "/")[end]
            table = read_mc_files(det, files)

            z = parse(Float64, split(p_name, "_")[5])
            φ = parse(Float64, split(p_name, "_")[7])

            h = append!(Histogram(binning), data)
            push!(histograms, HistogramWithPars(h, z = z, φ = φ))

            #
            push!(zs, z)
            push!(φs, φ)

        end
        outputs[det] =
            GeneralisedHistogram(histograms, z = _set2range(zs), φ = _set2range(φs))
    end
    return outputs
end
