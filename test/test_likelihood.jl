using CalibrationTemplateFits
using Test

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

end
