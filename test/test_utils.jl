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

    geds = Table(
        energy = VectorOfVectors([[100, 100, 200], [500]]),
        rawid = VectorOfVectors([[1, 2, 3], [1]]),
    )
    data = Table(geds = geds)
    @test CalibrationTemplateFits.get_data_histogram(1, data, 0:1000:1000) isa Histogram
    @test CalibrationTemplateFits.get_data_histogram(1, data, 0:1000:1000).weights == [2.0]
end
