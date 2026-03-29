using CalibrationTemplateFits
using Test
using StatsBase

@testset "test_likelihood" begin

    @test isapprox(
        CalibrationTemplateFits._poisson_likelihood([0], [2.0]),
        -2.0,
        atol = 1e-5,
    )
    @test isapprox(
        CalibrationTemplateFits._poisson_likelihood([0, 0], [2.0, 3.0]),
        -5.0,
        atol = 1e-5,
    )

    data = append!(Histogram(2600:100:2700), [])

    hists = HistogramWithPars[]

    h1 = append!(Histogram(2600:100:2700), [2610])
    h2 = append!(Histogram(2600:100:2700), [2620, 2610])

    push!(hists, HistogramWithPars(h1, par = 0))
    push!(hists, HistogramWithPars(h2, par = 1))

    model = GeneralisedHistogram(hists, par = 0:1:1)

    # likelihood should be just - model prediction
    @test isapprox(spectrum_likelihood(data, model, (par = 0.0, A = 1)), -1.0, atol = 1e-5)
    @test isapprox(spectrum_likelihood(data, model, (par = 1.0, A = 1)), -2.0, atol = 1e-5)
    @test isapprox(spectrum_likelihood(data, model, (par = 0.5, A = 1)), -1.5, atol = 1e-5)



end

@testset "test_build_likelihood" begin

    data = append!(Histogram(2600:100:2700), [])
    h1 = append!(Histogram(2600:100:2700), [2610])
    hists = [HistogramWithPars(h1, par = 0)]
    model = GeneralisedHistogram(hists, par = 0:1:0)

    data_hists = Dict(:det1 => data)
    models_dict = Dict(:det1 => model)

    lh = build_likelihood(data_hists, models_dict)
    @test lh !== nothing

    # Evaluate via the wrapped function (DensityInterface.LogFuncDensity stores log_f)
    @test isapprox(lh.log_f((par = 0.0, A = 1.0)), -1.0, atol = 1e-5)

    # livetime and n_sim scaling
    lh_scaled = build_likelihood(data_hists, models_dict, livetime = 2.0, n_sim = 2.0)
    @test isapprox(lh_scaled.log_f((par = 0.0, A = 1.0)), -1.0, atol = 1e-5)

end

@testset "test_build_prior" begin

    # Basic prior (no vary_fccd)
    prior = build_prior([:det1])
    @test prior !== nothing

    # Custom limits
    prior_lims =
        build_prior([:det1], zlims = (-5.0, 5.0), φlims = (-3.0, 3.0))
    @test prior_lims !== nothing

    # Prior with vary_fccd adds per-detector fccd parameters
    prior_fccd = build_prior([:det1, :det2], vary_fccd = true)
    @test prior_fccd !== nothing

end
