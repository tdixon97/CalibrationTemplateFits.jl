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


"""
    rebin_integer(h::Histogram, new_edges)

Rebin a 1D `StatsBase.Histogram` `h` to integer bin edges given by `new_edges`,
which can be an `AbstractRange` or a `Vector{Int}`.

Throws an error if `new_edges` are not integers, not sorted, or fall outside the
range of `h`.
"""
function rebin_integer(h::Histogram, new_edges)
    # Validate bin edges
    if !(new_edges isa AbstractVector{<:Integer} || new_edges isa AbstractRange{<:Integer})
        throw(ArgumentError("new_edges must be an AbstractRange or Vector of Ints"))
    end
    if !issorted(new_edges)
        throw(ArgumentError("new_edges must be sorted"))
    end

    # Check that new bin edges are within old edges
    old_edges = h.edges[1]
    if first(new_edges) < first(old_edges) || last(new_edges) > last(old_edges)
        throw(ArgumentError("new_edges must be within the range of the original histogram"))
    end

    # Create a new histogram and fill it
    rebinned = fit(Histogram, Float64[], new_edges)
    for (count, (lo, hi)) in zip(h.weights, zip(old_edges[1:(end-1)], old_edges[2:end]))
        mid = (lo+hi)/2
        bin_index = searchsortedlast(new_edges, mid)
        if (bin_index < length(new_edges)) & (bin_index > 0)
            rebinned.weights[bin_index] += count
        end
    end

    return rebinned
end
