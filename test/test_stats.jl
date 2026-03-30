
using CalibrationTemplateFits
using Test
using Distributions

@testset "test_normalised_poisson_residual" begin

    # Gaussian regime (mu >= 50): residual = (x - mu) / sqrt(mu)
    @test isapprox(normalised_poisson_residual(100.0, 100), 0.0, atol = 1e-10)
    @test isapprox(normalised_poisson_residual(100.0, 110), 1.0, atol = 1e-10)
    @test isapprox(normalised_poisson_residual(100.0, 90), -1.0, atol = 1e-10)

    # Low-mu regime: x == 0 and mode == 0 -> special-cased to return 0.0
    # (mu=0.3 -> mode=floor(0.3)=0, x=0 -> hits the x==0 && mode==0 branch)
    @test normalised_poisson_residual(0.3, 0) == 0.0
    # x < 1 is treated as x = 0; mode = floor(0.3) = 0 -> also hits x==0 && mode==0
    @test normalised_poisson_residual(0.3, 0.5) == 0.0

    # Low-mu regime: x < mode -> negative residual
    # mu=5, mode=5, x=1 -> sgn=-1
    @test normalised_poisson_residual(5.0, 1) < 0.0

    # Low-mu regime: x >= mode -> non-negative residual
    # mu=5, mode=5, x=10 -> sgn=1
    @test normalised_poisson_residual(5.0, 10) > 0.0

    # Exact values for known cases
    let μ = 5.0, x = 10
        prob = 1 - cdf(Poisson(μ), x - 1)
        expected = quantile(Normal(), 1 - prob)
        @test isapprox(normalised_poisson_residual(μ, x), expected, atol = 1e-10)
    end

    let μ = 5.0, x = 1
        prob = cdf(Poisson(μ), x)
        expected = -quantile(Normal(), 1 - prob)
        @test isapprox(normalised_poisson_residual(μ, x), expected, atol = 1e-10)
    end

    # Broadcasting over vectors
    mus = [100.0, 100.0, 5.0]
    obs = [110, 90, 0]
    res = normalised_poisson_residual(mus, obs)
    @test res isa AbstractVector
    @test length(res) == 3
    @test isapprox(res[1], 1.0, atol = 1e-10)
    @test isapprox(res[2], -1.0, atol = 1e-10)
    @test res[3] < 0.0   # mu=5, x=0 < mode=5 -> negative

    # Result type is Float64 for scalar inputs
    @test normalised_poisson_residual(100.0, 100) isa Float64

end
