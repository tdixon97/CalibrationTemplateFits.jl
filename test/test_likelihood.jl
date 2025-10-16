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
    @test isapprox(spectrum_likelihood(data, model, (par = 0.0,)), -1.0, atol = 1e-5)
    @test isapprox(spectrum_likelihood(data, model, (par = 1.0,)), -2.0, atol = 1e-5)
    @test isapprox(spectrum_likelihood(data, model, (par = 0.5,)), -1.5, atol = 1e-5)



end
