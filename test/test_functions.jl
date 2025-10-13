
using Test

@testset "test_activeness" begin

    @test piecewise_linear_activeness(2, fccd = 1, dlf = 0.5) == 1.0
    @test piecewise_linear_activeness(0.2, fccd = 1, dlf = 0.5) == 0.0
    @test piecewise_linear_activeness(0.75, fccd = 1, dlf = 0.5) == 0.5

end
