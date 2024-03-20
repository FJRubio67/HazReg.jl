# Function to standardise a vector.
function standardise(x)
    out = (x .- mean(x)) ./ std(x)
     return out
 end

# Standard definitions for hazard rates that work for every distributions at once: 
loghaz(d::UnivariateDistribution,t::Real) = Distributions.logpdf(d,t) - Distributions.logccdf(d,t)
haz(d::UnivariateDistribution,t::Real) = exp(loghaz(d,t))
cumhaz(d::UnivariateDistribution,t::Real) = -Distributions.logccdf(d,t)