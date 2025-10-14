
using LegendHDF5IO
using StatsBase
using ProgressMeter
using TypedTables
using StatsBase

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
    get_data_histogram(det::Symbol,data::Table,range::E)

Extract the histogram of the data.
"""
function get_data_histogram(det::Symbol, data::Table, range)
    r = chmap[det].daq.rawid

    # flatten and get energy and rawid
    rawid = flatview(data.geds.rawid)
    energy = flatview(data.geds.energy)
    energy_sel = energy[rawid .== r]

    return append!(Histogram(range), energy_sel)
end

"""Convert a set to a range"""
function _set2range(set::Set)

    vector = sort(collect(set))

    range = vector[1]:diff(vector[1]):vector[end]

    if collect(range)!=vector
        error("vector is not evenly spaced, cannot convert to StepRangeLen")
    end
    range
end

"""
    read_mc_models(mc_path::String, label::String,binning::Union{AbstractVector,AbstractRange})

Read the MC files from `mc_path`,saving the results to the `GeneralisedHistogram`` models.
"""
function read_mc_models(
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
        phis = Set([])

        for folder in positions

            files = glob("*", folder)
            p_name = split(folder, "/")[end]
            table = read(det, files)

            z = parse(Float64, split(p_name, "_")[5])
            phi = parse(Float64, split(p_name, "_")[7])

            h = append!(Histogram(binning), data)
            push!(histograms, HistogramWithPars(h, z = z, phi = phi))

            #
            push!(zs, z)
            push!(phis, phi)

        end
        outputs[det] =
            GeneralisedHistogram(histograms, z = _set2range(zs), phi = _set2range(phis))
    end
    return outputs
end


"""
    create_generalised_hist(mc::Table,binning::AbstractRange; vary_fccd::Bool = false, vary_dlf::Bool = false)

Build the generalised histogram from the MC.
"""
function create_generalised_hist(
    mc::Vector,
    binning,
    z_range::AbstractRange,
    phi_range::AbstractRange,
)
    hs = HistogramWithPars[]
    for z in collect(z_range)
        for phi in collect(phi_range)

            data = nothing
            for m in mc
                if (m.z == z) & (m.phi == phi)
                    data = m.data
                end
            end
            hist = HistogramWithPars(append!(Histogram(range), data), z = z, phi = phi)
            push!(hs, hist)
        end
    end
    ghist = GeneralisedHistogram(hs, z = z_range, phi = phi_range)
    return ghist

end

"""
    parse_binning(s::AbstractString) -> Edges

Parses a string describing either a range (`"a:b:c"`) or a vector (`"[x₁, x₂, ...]"`)
and wraps the result in a `Histograms.Edges` object.

"""
function parse_binning(s::AbstractString)
    s = strip(s)
    if startswith(s, "[") && endswith(s, "]")
        nums = split(s[2:(end-1)], ',')
        edges = parse.(Float64, strip.(nums))
        return Edges(edges)
    else
        parts = split(s, ':')
        if length(parts) != 3
            error("Expected range of form a:b:c, got $s")
        end
        a, b, c = parse.(Float64, parts)
        return Edges(a:b:c)
    end
end
