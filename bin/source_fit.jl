

using Pkg
Pkg.activate(".") # activate the environment
Pkg.instantiate()

using CalibrationTemplateFits
using BAT


function main()

    @info "running using ", Base.Threads.nthreads(), " threads"

    dets = get_detectors()

    mc = read_mc(input_mc)

    # get the models
    models = GeneralisedHistogram[]

    for det in dets
        push!(models, get_model(mc, det; model_kwargs...))
    end

    # get the data hist
    data_table = read_data(input)

    data_histograms = get_data_hist(data_table; model_kwargs...)
    # build the likelihood
    posterior = build_posterior(data, models; priors_kwargs...)

    # sample
    samples =
        bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(); fit_kwargs...)).result

    # save
    save_result(samples, output)
end
