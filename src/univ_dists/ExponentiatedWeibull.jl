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
        sigma,nu,gamma = promote(eltype.((sigma,nu,gamma))...)
        new{T}(T(sigma), T(nu), T(gamma), Weibull(nu,sigma))
    end
end
params(d::ExponentiatedWeibull) = (d.sigma,d.nu,d.gamma) # should be such that D(params(d)) == d

@distr_support ExponentiatedWeibull 0.0 Inf

function rand(rng::AbstractRNG, d::ExponentiatedWeibull)
    p = rand(rng) ^ (1 / d.gamma)
    return quantile(d.W, p)
end
function logpdf(d::ExponentiatedWeibull, t::Real)
    return log(d.gamma) + (d.gamma - 1) * logcdf(d.W, t) + logpdf(d.W, t)
end
# pdf is supplied by default from the interface. 
function logcdf(d::ExponentiatedWeibull, t::Real)
    return d.gamma * logcdf(d.W, t)
end
cdf(d::MyDistribution, t::Real) = exp(logcdf(d,t))
function quantile(d::ExponentiatedWeibull, p::Real)
    return quantile(d.W, p .^ (1 / d.gamma))
end