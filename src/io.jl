
using LegendHDF5IO
using StatsBase
using ProgressMeter
using TypedTables

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

"""
    read_mc(mc_path::String, label::String)
"""
function read_mc(mc_path::String, label::String = "sis-1_source-2")
    positions = glob("$(label)_$pos*", mc_path)

    n_prim = nothing
    outputs = Dict{Symbol,Any}(ch => [] for ch in keys(chmap) if chmap[ch][:system]=="geds")

    @showprogress for det in keys(chmap)

        if (chmap[det][:system]!="geds")
            continue
        end
        for folder in positions

            files = glob("*", folder)
            p_name = split(folder, "/")[end]
            table = read(det, files)

            z = parse(Float64, split(p_name, "_")[5])
            phi = parse(Float64, split(p_name, "_")[7])

            push!(outputs[det], (data = table, phi = phi, z = z))

        end
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
