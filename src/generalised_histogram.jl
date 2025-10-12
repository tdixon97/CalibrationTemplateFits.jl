using StatsBase
using Interpolations
using Base.Iterators: product

"""A histogram with associated parameters.

Parameters
* hist: the `Histogram` itself
* pars: a set of parameters associated with the histogram.
"""
struct HistogramWithPars{T<:Real,N}
    hist::Histogram{T,N}
    pars::NamedTuple
end

"""
    HistogramWithPars(hist::Histogram{T,N}; kwargs...)

Convienience constructor for the `HistogramWithPars` from keyword arguments.
"""
function HistogramWithPars(hist::Histogram{T,N}; kwargs...) where {T,N}
    return HistogramWithPars{T,N}(hist, NamedTuple(kwargs))
end

"""
    _find_histogram(histograms::AbstractVector,pars::NamedTuple)

Find the histogram corresponding to a given set of parameters.
"""
function _find_histogram(histograms::AbstractVector{<:HistogramWithPars}, pars::NamedTuple)

    for h in histograms

        if h.pars == pars
            return h
        end
    end

    throw(ArgumentError("Histogram with parameter $pars not found."))

end

"""
    get_counts_grid(idx::Int,histograms::AbstractVector,grid::NamedTuple)

Extract the counts evaluated on the `grid`. This is done by iterating over all the combinations
of fields in the `grid`, and for each finding the corresponding histogram. We then extract
the bin content corresponding to `idx`.
"""
function _get_counts_grid(
    idx::Int,
    histograms::AbstractVector{<:HistogramWithPars},
    grid::G,
) where {G<:NamedTuple}

    ntemp = [length(gt) for gt in values(grid)]
    n = prod(ntemp)

    counts = Vector{Float64}(undef, n)

    for (i, vals) in enumerate(Iterators.product(grid...))
        named = NamedTuple(keys(grid) .=> vals)

        # find the histogram
        h_tmp = _find_histogram(histograms, named)

        # save to the flat array
        counts[i] = h_tmp.hist.weights[idx]
    end

    # reshape into a grid
    return reshape(counts, ntemp...)

end

raw"""The `GeneralisedHistogram` type represents
a binned model where the counts in every bin
is a conditional distribution. Bin-by-bin
interpolation is performed to extract the
weights.

Currently only one dimensional histograms with uniformly
spaced grids are supported.

Mathematically for every bin $i$, the `GeneralisedHistogram`
stores the conditional distribution

```math
    p(N_i|\theta)
```
where $N_i$ is the bin content for $i$ and $\theta$ are parameters.

This can be used to construct generalised binned models where
the model predictions are functions of nuisance parameters.

# Fields
* interpolators: The interpolator object per bin
* edges: The bin edges
* grid: A dictonary of the grid for nuisance parameters

# Examples

```julia-repl
julia> histograms = [HistogramWithPars(Histogram(0:1:4000),fccd=1.),
                    HistogramWithPars(Histogram(0:1:4000),fccd=2.)]
julia> ghist = GeneralisedHist(histograms,fccd = 1.:1:2.)
julia> get_bin_content(ghist,fccd = 1.5)
0.
```

An arbitrary number of parameters is supported!
"""
struct GeneralisedHistogram
    interpolators::Vector{Any}
    edges::Any
    grid::NamedTuple
end

"""
    GeneralisedHistogram(histograms::AbstractVector; kwargs...)

Construct the `GeneralisedHistogram` using keyword arguments.
"""
function GeneralisedHistogram(histograms::AbstractVector; kwargs...)
    grid = NamedTuple(kwargs)
    return GeneralisedHistogram(histograms, grid)
end

"""
    function GeneralisedHistogram(histograms::AbstractVector, grid::G)

Construct the generalised histogram from a list of `HistogramWithPars` and an associated parameter grid.

`histograms` should be a list of `HistogramWithPars` for each parameter value in the `grid`, which is a `NamedTuple`
of the grid spacing for every parameter.

"""
function GeneralisedHistogram(histograms::AbstractVector, grid::G) where {G}
    edges = histograms[1].hist.edges[1]
    size = length(histograms[1].hist.weights)

    # check the edges
    for hist in histograms
        if edges != hist.hist.edges[1]
            throw(ArgumentError("All histograms must have the same edges"))
        end
    end

    # build interpolators
    interpolators = Vector{Any}()
    for i = 1:size
        counts_grid = _get_counts_grid(i, histograms, grid)
        push!(interpolators, interpolate(counts_grid, BSpline(Linear())))
    end

    return GeneralisedHistogram(interpolators, edges, grid)
end

"""
    grid_value(range::AbstractRange,point::Real)->Real

Extract the value of the `point` in the integer coordinates of the
`range`. These coordinates are "1" for the first point and `length(step)`
for the last.
"""
function grid_value(range::AbstractRange, point::Real)

    step_size = step(range)

    return ((point - first(range))/step_size)+1
end

"""
    get_normalised_par_values(grid::NamedTuple,pars::NamedTuple)

Extract the parameter values on the grid.
"""
function get_normalised_par_values(grid::NamedTuple, pars::NamedTuple)

    norm_param_values = Tuple(grid_value(grid[k], pars[k]) for k in keys(grid))
    if (length(norm_param_values)==1)
        return norm_param_values[1]
    else
        return norm_param_values
    end
end

"""
    get_weights(hist::GeneralisedHistogram; kwargs...) -> Vector

Get the weights evaluating the generalising histogram at the kwargs (nuisance parameters).
"""
function get_weights(hist::GeneralisedHistogram; kwargs...)
    parameters = NamedTuple(kwargs)

    nbins = length(hist.interpolators)
    weights = zeros(nbins)

    # normalise the parameter value onto the unit grid
    norm_param_values = get_normalised_par_values(hist.grid, parameters)

    for i = 1:nbins
        weights[i] = hist.interpolators[i](norm_param_values...)
    end

    return weights
end

"""
    get_histogram(hist::GeneralisedHistogram; kwargs...) -> Histogram
"""
function get_histogram(hist::GeneralisedHistogram; kwargs...)

    histogram = fit(Histogram{Float64}, Float64[], hist.edges)
    histogram.weights .= Float64.(get_weights(hist; kwargs...))

    return histogram
end


"""
    get_bin_content(idx::Int,hist::GeneralisedHistogram; kwargs)-> T

Get the bin content for bin `idx` evaluating the generalising histogram.
"""
function get_bin_content(idx::Int, hist::GeneralisedHistogram; kwargs...)

    parameters = NamedTuple(kwargs)

    # normalise the parameter value onto the unit grid
    # TODO: test the order
    norm_param_values = get_normalised_par_values(hist.grid, parameters)

    return hist.interpolators[idx](norm_param_values...)
end
