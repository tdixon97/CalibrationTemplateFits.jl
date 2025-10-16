
using StatsBase
using Distributions
using BAT
using DensityInterface

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
    params::NamedTuple,
)
    # get model predictions
    mu = forward_model(model, params)
    _poisson_likelihood(data.weights, mu)

end

"""
    build_likelihood(data_hists::Dict, models::Dict)

Build the likelihood summing over all the spectra in the
`data_hists` and models defined by `models`. In each case
a Poisson likelihood is employed.
"""
function build_likelihood(data_hists::Dict, models::Dict)

    logfuncdensity(
        function (params::NamedTuple)

            ll_value = 0

            # loop over spectra
            for ch in keys(data_hists)
                ll_value += spectrum_likelihood(data_hists[ch], models[ch], params)
            end

            return ll_value
        end,
    )
end
