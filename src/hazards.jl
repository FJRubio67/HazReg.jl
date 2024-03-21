#= 
=========================================================================
Hazard Functions, PDFs, CDFs
========================================================================= 
=#

#= LogNormal hazard function.
t: positive argument

mu: log-mean parameter

sigma: scale parameter

logh: log scale (true or false)
=#

function hLogNormal(t, mu, sigma, logh::Bool=false)
    lpdf0 = logpdf.(LogNormal(mu, sigma), t)
    ls0 = logccdf.(LogNormal(mu, sigma), t)
    val = lpdf0 .- ls0
    if logh
        return val
    else
        return exp.(val)
    end
end

#= LogNormal cumulative hazard function.
t: positive argument

mu: log-mean parameter

sigma: scale parameter

=#

function chLogNormal(t, mu, sigma)
    H0 = -logccdf.(LogNormal(mu, sigma), t)
    return H0
end


#= LogLogistic hazard function.
t: positive argument

mu: log-mean parameter

sigma: scale parameter

logh: log scale (true or false)
=#

function hLogLogistic(t, mu, sigma, logh::Bool=false)
    lt = log.(t)
    lpdf0 = logpdf.(Logistic(mu, sigma), lt) .- lt
    ls0 = logccdf.(Logistic(mu, sigma), lt)
    val = lpdf0 .- ls0
    if logh
        return val
    else
        return exp.(val)
    end
end

#= LogLogistic cumulative hazard function.
t: positive argument

mu: log-mean parameter

sigma: scale parameter
=#

function chLogLogistic(t, mu, sigma)
    lt = log.(t)
    H0 = -logccdf.(Logistic(mu, sigma), lt)
    return H0
end


#= Weibull hazard function.
t: positive argument

shape: shape parameter

scale: scale parameter

logh: log scale (true or false)
=#

function hWeibull(t, shape, scale, logh::Bool=false)
    lpdf0 = logpdf.(Weibull(shape, scale), t)
    ls0 = logccdf.(Weibull(shape, scale), t)
    val = lpdf0 .- ls0
    if logh
        return val
    else
        return exp.(val)
    end
end

#= Weibull cumulative hazard function.
t: positive argument

shape: shape parameter

scale: scale parameter
=#

function chWeibull(t, shape, scale)
    H0 = -logccdf.(Weibull(shape, scale), t)
    return H0
end

#= Gamma hazard function.
t: positive argument

shape: shape parameter

scale: scale parameter

logh: log scale (true or false)
=#

function hGamma(t, shape, scale, logh::Bool=false)
    lpdf0 = logpdf.(Gamma(shape, scale), t)
    ls0 = logccdf.(Gamma(shape, scale), t)
    val = lpdf0 .- ls0
    if logh
        return val
    else
        return exp.(val)
    end
end

#= Gamma cumulative hazard function.
t: positive argument

shape: shape parameter

scale: scale parameter
=#

function chGamma(t, shape, scale)
    H0 = -logccdf.(Gamma(shape, scale), t)
    return H0
end