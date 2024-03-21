"""
    LogLogistic(mu,sigma)

To be described... 


```math
f(x; parameters) = ...
```

"""
struct LogLogistic{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    mu::T
    sigma::T
    function LogLogistic(mu,sigma)
        T = promote_type(Float64, eltype.((mu,sigma))...)
        return new{T}(T(mu), T(sigma))
    end
end
LogLogistic() = LogLogistic(1,1)
Distributions.params(d::LogLogistic) = (d.mu,d.sigma)
Distributions.@distr_support LogLogistic 0.0 Inf
function loghaz(d::LogLogistic, t::Real)
    lt = log.(t)
    lpdf0 = Distributions.logpdf.(Logistic(d.mu, d.sigma), lt) .- lt
    ls0 = Distributions.logccdf.(Logistic(d.mu, d.sigma), lt)
    return lpdf0 .- ls0
end
function cumhaz(d::LogLogistic,t::Real)
    lt = log.(t)
    return -Distributions.logccdf.(Logistic(d.mu, d.sigma), lt)
end

# We could define more : 

# function Distributions.rand(rng::AbstractRNG, d::LogLogistic)

# end
# function Distributions.logpdf(d::LogLogistic, t::Real)
    
# end
# function Distributions.logccdf(d::LogLogistic, t::Real)
    
# end
# Distributions.ccdf(d::LogLogistic, t::Real) = exp(Distributions.logccdf(d,t))
# Distributions.cdf(d::LogLogistic, t::Real) = 1 - Distributions.ccdf(d,t)
# function Distributions.quantile(d::LogLogistic, p::Real)
    
# end
