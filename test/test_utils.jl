using CalibrationTemplateFits
using Test
using TypedTables
using ArraysOfArrays
using StatsBase

@testset "test_set2range" begin

    @test CalibrationTemplateFits._set2range(Set([1, 2, 3.0]))==1.0:1:3.0
end

@testset "test_binning" begin

    @test parse_binning("0.:10:20") == 0:10:20
    @test parse_binning("[1,2,3]") == [1, 2, 3]
end


@testset "test_histograms" begin

    quality = Table(is_good_channel = VectorOfVectors([[true, true, true], [true]]))
    geds = Table(
        energy = VectorOfVectors([[100, 100, 200], [500]]),
        rawid = VectorOfVectors([[1, 2, 3], [1]]),
        quality = quality,
    )

    data = Table(geds = geds)

    @test CalibrationTemplateFits.get_data_histogram(1, data, 0:1000:1000) isa Histogram
    @test CalibrationTemplateFits.get_data_histogram(1, data, 0:1000:1000).weights == [2.0]
end



@testset "rebin_integer" begin
    h = fit(Histogram, rand(0:9, 1000), 0:1:10)

    # 1. Basic rebinning
    @testset "Basic rebinning" begin
        new_edges = 0:2:10
        h2 = CalibrationTemplateFits.rebin_integer(h, new_edges)
        @test h2.edges[1] == collect(new_edges)
        @test sum(h2.weights) ≈ sum(h.weights)
        @test length(h2.weights) == length(new_edges) - 1
    end

    # 2. Identity (same edges)
    @testset "Identity" begin
        h2 = CalibrationTemplateFits.rebin_integer(h, 0:1:10)
        @test all(isapprox.(h.weights, h2.weights))
    end


    # 3. One-bin edge case
    @testset "One-bin edge case" begin
        h2 = CalibrationTemplateFits.rebin_integer(h, 0:10:10)
        @test sum(h2.weights) ≈ sum(h.weights)
        @test length(h2.weights) == 1
    end
end


@testset "extract_grid_values" begin
    path = joinpath(@__DIR__, "test_files")

    models = read_models(
        [:det1, :det2, :det3],
        [
            path*"/hit_files/z_-1_phi_0",
            path*"/hit_files/z_1_phi_0",
            path*"/hit_files/z_-1_phi_1",
            path*"/hit_files/z_1_phi_1",
        ],
        2600:20:2620,
        r".*z_([-\d.]+)_phi_([-\d.]+)",
    )

    @test CalibrationTemplateFits.extract_grid_values(models, :det1) ==
          ((-1.0, 1.0), (0.0, 1.0))
    @test CalibrationTemplateFits.extract_grid_values(models, :det2) ==
          ((-1.0, 1.0), (0.0, 1.0))

end
