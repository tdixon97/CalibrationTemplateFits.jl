using Distributions

"""
    normalised_poisson_residual(mu, obs)

Compute normalised Poisson residuals.

For `mu > 50`, a Gaussian approximation is used.
When `mu < 50`, the Poisson tail probability is calculated
from the mode, used to derive a p-value, and converted into a residual.

# Arguments
- `mu`: model predictions (scalar or array)
- `obs`: observed data (scalar or array)

# Returns
A scalar or array of normalized residuals with the same shape as `mu` and `obs`.
"""
function normalised_poisson_residual(mu, obs)
    broadcast(mu, obs) do μ, x
        x < 1 && (x = 0)
        if μ < 50
            mode = floor(Int, μ)
            if x == 0 && mode == 0
                return 0.0
            end

            if x < mode
                prob = cdf(Poisson(μ), x)
                sgn = -1
            else
                prob = 1 - cdf(Poisson(μ), x - 1)
                sgn = 1
            end
            return sgn * quantile(Normal(), 1 - prob)
        else
            return (x - μ) / sqrt(μ)
        end
    end
end
