using CalibrationTemplateFits
using Test

@testset "test_set2range" begin

    @test CalibrationTemplateFits._set2range(Set([1, 2, 3.0]))==1.0:1:3.0
end

@testset "test_binning" begin

    @test parse_binning("0.:10:20") == 0:10:20
    @test parse_binning("[1,2,3]") == [1, 2, 3]
end
