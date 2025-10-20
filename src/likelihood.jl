
using StatsBase
using Distributions
using BAT
using DensityInterface
using IntervalSets

"""
    poisson_likelihood(obs_weight::Vector, pred_weight::Vector)

Get the likelihood for a series of independent Poisson distributed random
variables.
"""
function _poisson_likelihood(obs_weight::AbstractVector, pred_weight::AbstractVector)

    return sum(logpdf.(Poisson.(pred_weight .+ eps.(pred_weight)), obs_weight))
end


"""
    spectrum_likelihood(data::Histogram,model::GeneralisedHistogram,params::NamedTuple)

Get the likelihood for a given spectrum, based on the `data` histogram and the `model` which is
in this case a `GeneralisedHistogram`.
"""
function spectrum_likelihood(
    data::Histogram,
    model::GeneralisedHistogram,
    params::NamedTuple;
    norm::Real = 1,
)
    # get model predictions
    mu = get_weights(model; params...)*params.A*norm
    _poisson_likelihood(data.weights, mu)

end

"""
    build_likelihood(data_hists::Dict, models::Dict)

Build the likelihood summing over all the spectra in the
`data_hists` and models defined by `models`. In each case
a Poisson likelihood is employed.
"""
function build_likelihood(
    data_hists::Dict,
    models::Dict;
    livetime::Real = 1,
    n_sim::Real = 1,
)

    logfuncdensity(
        function (params::NamedTuple)

            ll_value = 0

            # loop over spectra
            for ch in keys(data_hists)
                ll_value += spectrum_likelihood(
                    data_hists[ch],
                    models[ch],
                    params,
                    norm = livetime/n_sim,
                )
            end

            return ll_value
        end,
    )
end

function build_prior(dets::AbstractVector; vary_fccd::Bool = false)

    dist = if !vary_fccd
        distprod(A = 0.0 .. 3000.0, z = -20.0 .. 20.0, φ = -6.0 .. 6.0)
    else
        dists = Dict(:A=>0.0 .. 3000, :z=>-80 .. -40.0, :φ => -6.0 .. 6.0)
        for det in dets
            dists[Symbol("fccd_$det")] = 0 .. 2
        end
        distprod(; dists...)
    end
    dist
end
