"""
    ExponentiatedWeibull(sigma,nu,gamma)

The *ExponentiatedWeibull distribution* with scale `sigma`, shape `nu` and second shape `gamma` has probability density function 

```math
f(x; parameters) = ...
```

More details and examples of usage could be provided in this docstring.

Maybe this distribution could simply be constructed from a transformation of the original Weibull ? 

References: 
* [Link to my reference so that people understand what it is](https://myref.com) 
"""
struct ExponentiatedWeibull{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    sigma::T
    nu::T
    gamma::T
    W::Weibull{T} # underlying Weibull distribution. 
    function ExponentiatedWeibull(sigma,nu,gamma)
        T = promote_type(Float64, eltype.((sigma,nu,gamma))...)
        return new{T}(T(sigma), T(nu), T(gamma), Weibull(nu,sigma))
    end
end
ExponentiatedWeibull() = ExponentiatedWeibull(1,1,1)
Distributions.params(d::ExponentiatedWeibull) = (d.sigma,d.nu,d.gamma) # should be such that D(params(d)) == d

Distributions.@distr_support ExponentiatedWeibull 0.0 Inf

function Distributions.rand(rng::AbstractRNG, d::ExponentiatedWeibull)
    p = Base.rand(rng) ^ (1 / d.gamma)
    return Distributions.quantile(d.W, p)
end
function Distributions.logpdf(d::ExponentiatedWeibull, t::Real)
    return log(d.gamma) + (d.gamma - 1) * Distributions.logcdf(d.W, t) + Distributions.logpdf(d.W, t)
end
# pdf is supplied by default from the interface. 
function Distributions.logcdf(d::ExponentiatedWeibull, t::Real)
    return d.gamma * Distributions.logcdf(d.W, t)
end
Distributions.cdf(d::ExponentiatedWeibull, t::Real) = exp(Distributions.logcdf(d,t))
function Distributions.quantile(d::ExponentiatedWeibull, p::Real)
    return Distributions.quantile(d.W, p .^ (1 / d.gamma))
end