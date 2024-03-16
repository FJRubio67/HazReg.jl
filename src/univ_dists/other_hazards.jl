
# Standard definitions for hazard rates that work for every distributions at once: 
loghaz(d::UnivariateDistribution,t::Real) = logpdf(d,t) - logccdf(d,t)
haz(d::UnivariateDistribution,t::Real) = exp(loghaz(d,t))
cumhaz(d::UnivariateDistribution,t::Real) = -logccdf(d,t)
