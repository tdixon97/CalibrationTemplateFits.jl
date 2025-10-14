
using LegendHDF5IO
using StatsBase
using ProgressMeter
using TypedTables
using StatsBase
using CalibrationTemplateFits

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
    read_models(mc_path::String, label::String,binning::Union{AbstractVector,AbstractRange})

Read the MC files from `mc_path`, saving the results to the `GeneralisedHistogram`` models.
"""
function read_models(
    mc_path::String,
    label::String,
    binning::Union{AbstractVector,AbstractRange},
)
    positions = glob("$(label)_$pos*", mc_path)

    n_prim = nothing
    outputs = Dict()

    @showprogress for det in keys(chmap)

        if (chmap[det][:system]!="geds")
            continue
        end
        histograms = HistogramWithPars[]

        zs = Set([])
        φs = Set([])

        for folder in positions

            files = glob("*", folder)
            p_name = split(folder, "/")[end]
            table = read(det, files)

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
