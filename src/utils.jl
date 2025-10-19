using TypedTables
using ArraysOfArrays

"""Convert a set to a range"""
function _set2range(set::Set)

    vector = Vector(sort(collect(set)))

    range = vector[1]:diff(vector)[1]:vector[end]

    if collect(range)!=vector
        error("vector is not evenly spaced, cannot convert to StepRangeLen")
    end
    range
end

"""
    parse_binning(s::AbstractString)

Parses a string describing either a range (`"a:b:c"`) or a vector (`"[x₁, x₂, ...]"`)
and wraps the result in a `Histograms.Edges` object.

"""
function parse_binning(s::AbstractString)
    s = strip(s)
    if startswith(s, "[") && endswith(s, "]")
        nums = split(s[2:(end-1)], ',')
        edges = parse.(Float64, strip.(nums))
        return edges
    else
        parts = split(s, ':')
        if length(parts) != 3
            error("Expected range of form a:b:c, got $s")
        end
        a, b, c = parse.(Float64, parts)
        return a:b:c
    end
end


"""
    get_data_histogram(det::Symbol,data::Table,range)

Extract the histogram of the data.
"""
function get_data_histogram(rawid::Int, data::Table, range)

    # flatten and get energy and rawid
    r = flatview(data.geds.rawid)
    energy = flatview(data.geds.energy)
    is_good = flatview(data.geds.quality.is_good_channel)

    energy_sel = energy[(r .== rawid) .& is_good]

    return append!(Histogram(range), energy_sel)
end
