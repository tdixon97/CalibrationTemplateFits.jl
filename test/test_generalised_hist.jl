using CalibrationTemplateFits
using Test
using Random
using Distributions
using StatsBase

@testset "test_basic" begin

    # test converting to grid coordinates
    @test CalibrationTemplateFits.grid_value(0.0:0.1:5, 0.0) == 1.0
    @test CalibrationTemplateFits.grid_value(0.0:0.1:5, 5.0) == 51.0
    @test CalibrationTemplateFits.grid_value(-2.5:0.1:2.5, 0) == 26.0

end

@testset "test_generalised_hist" begin

    # construct the inputs
    hists = HistogramWithPars[]
    for reso = 0.1:1.0:10
        data = rand(Normal(1000, reso), 10000)
        hist_tmp = HistogramWithPars(append!(Histogram(0:1.0:4000), data), reso = reso) # note the resolution is stored
        push!(hists, hist_tmp)
    end
    ghist = GeneralisedHistogram(hists, reso = 0.1:1.0:10)

    @test hists[1] isa HistogramWithPars
    @test ghist isa GeneralisedHistogram

    # test it throws with different histos
    hists_bad = HistogramWithPars[]
    push!(hists_bad, HistogramWithPars(append!(Histogram(0:1.0:4000), [1, 2, 3]), reso = 1))
    push!(hists_bad, HistogramWithPars(append!(Histogram(0:2.0:4000), [1, 2, 3]), reso = 2))


    @test_throws ArgumentError GeneralisedHistogram(hists_bad, reso = 1:1:2)

    # if we try to construct with the wrong range we get an error
    @test_throws ArgumentError GeneralisedHistogram(hists, sigma = 0.1:0.15:2)

    @testset "test_methods" begin

        # bins far from the peak should be 0
        @test CalibrationTemplateFits.get_bin_content(5, ghist, reso = 0.5)==0.0
        @test CalibrationTemplateFits.get_weights(ghist, reso = 0.5) isa Vector
        @test CalibrationTemplateFits.get_histogram(ghist, reso = 0.5) isa Histogram

        # everything in range should work
        for val = 0.1:0.05:10
            @test CalibrationTemplateFits.get_bin_content(5, ghist, reso = 0.5) == 0.0

        end

        # anything outside should break
        @test_throws BoundsError CalibrationTemplateFits.get_bin_content(
            5,
            ghist,
            reso = -5,
        )


    end
end

@testset "test_find_histogram" begin

    hists = [
        HistogramWithPars(append!(Histogram(0:1.0:10), Float64[]), par = 0.0),
        HistogramWithPars(append!(Histogram(0:1.0:10), Float64[]), par = 1.0),
    ]

    # Should find the correct histogram
    h = CalibrationTemplateFits._find_histogram(hists, (par = 0.0,))
    @test h.pars == (par = 0.0,)

    h = CalibrationTemplateFits._find_histogram(hists, (par = 1.0,))
    @test h.pars == (par = 1.0,)

    # Should throw when histogram not found
    @test_throws ArgumentError CalibrationTemplateFits._find_histogram(hists, (par = 2.0,))

end

@testset "test_get_normalised_par_values" begin

    # 1D grid
    grid_1d = (par = 0.0:1.0:4.0,)

    @test CalibrationTemplateFits.get_normalised_par_values(
        grid_1d,
        (par = 0.0,),
        Val(1),
    ) == (1.0,)
    @test CalibrationTemplateFits.get_normalised_par_values(
        grid_1d,
        (par = 2.0,),
        Val(1),
    ) == (3.0,)

    # 2D grid: returns a Tuple
    grid_2d = (z = -1.0:1.0:1.0, φ = 0.0:1.0:1.0)
    result = CalibrationTemplateFits.get_normalised_par_values(
        grid_2d,
        (z = 0.0, φ = 0.5),
        Val(2),
    )
    @test result isa Tuple
    # grid_value(range, point) = (point - first(range)) / step(range) + 1
    # z: (0.0 - (-1.0)) / step(-1.0:1.0:1.0) + 1 = 1.0/1.0 + 1 = 2.0
    @test isapprox(result[1], 2.0, atol = 1e-10)
    # φ: (0.5 - 0.0) / step(0.0:1.0:1.0) + 1 = 0.5/1.0 + 1 = 1.5
    @test isapprox(result[2], 1.5, atol = 1e-10)

    # test a different order

    result = CalibrationTemplateFits.get_normalised_par_values(
        grid_2d,
        (b = 10, z = 0.0, φ = 0.5),
        Val(2),
    )

    @test result isa Tuple
    # grid_value(range, point) = (point - first(range)) / step(range) + 1
    # z: (0.0 - (-1.0)) / step(-1.0:1.0:1.0) + 1 = 1.0/1.0 + 1 = 2.0
    @test isapprox(result[1], 2.0, atol = 1e-10)
    # φ: (0.5 - 0.0) / step(0.0:1.0:1.0) + 1 = 0.5/1.0 + 1 = 1.5
    @test isapprox(result[2], 1.5, atol = 1e-10)

end

@testset "test_2d_generalised_hist" begin

    # Build a 2x2 grid of histograms (z x φ)
    hists_2d = HistogramWithPars[]
    for z in [-1.0, 1.0], φ in [0.0, 1.0]
        h = append!(Histogram(0:1.0:10), Float64[])
        # Set the first bin to a known value based on parameters
        h.weights[1] = z + φ + 2.0  # so we can verify interpolation
        push!(hists_2d, HistogramWithPars(h, z = z, φ = φ))
    end

    ghist_2d = GeneralisedHistogram(hists_2d, z = -1.0:2.0:1.0, φ = 0.0:1.0:1.0)
    @test ghist_2d isa GeneralisedHistogram

    # At grid point (-1, 0): weight = -1 + 0 + 2 = 1.0
    @test isapprox(
        CalibrationTemplateFits.get_bin_content(1, ghist_2d, z = -1.0, φ = 0.0),
        1.0,
        atol = 1e-10,
    )

    # At grid point (1, 1): weight = 1 + 1 + 2 = 4.0
    @test isapprox(
        CalibrationTemplateFits.get_bin_content(1, ghist_2d, z = 1.0, φ = 1.0),
        4.0,
        atol = 1e-10,
    )

    # Interpolated at center (0, 0.5): weight should be average = 2.5
    @test isapprox(
        CalibrationTemplateFits.get_bin_content(1, ghist_2d, z = 0.0, φ = 0.5),
        2.5,
        atol = 1e-10,
    )

    @test CalibrationTemplateFits.get_weights(ghist_2d, z = 0.0, φ = 0.5) isa Vector

end
