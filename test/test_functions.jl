
using Test

@testset "test_activeness" begin

    @test piecewise_linear_activeness(2, fccd = 1, dlf = 0.5) == 1.0
    @test piecewise_linear_activeness(0.2, fccd = 1, dlf = 0.5) == 0.0
    @test piecewise_linear_activeness(0.75, fccd = 1, dlf = 0.5) == 0.5

    # Boundary: d == fccd falls into the else branch -> returns 1.0
    @test piecewise_linear_activeness(1.0, fccd = 1.0, dlf = 0.5) == 1.0

    # Boundary: d == dlf * fccd enters the else branch;
    # dl = fccd*dlf = 0.5, so (d - dl)/(fccd - dl) = (0.5 - 0.5)/(1.0 - 0.5) = 0.0
    @test piecewise_linear_activeness(0.5, fccd = 1.0, dlf = 0.5) == 0.0

    # Quarter-way through the transition region
    @test isapprox(
        piecewise_linear_activeness(0.625, fccd = 1.0, dlf = 0.5),
        0.25,
        atol = 1e-10,
    )

end
