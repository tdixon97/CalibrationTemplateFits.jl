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
