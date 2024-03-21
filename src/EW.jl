#=
Exponentiated Weibull Distribution: Hazard, cumulative hazard,
probability density function, random number generation, and survival function
=#

#= Exponentiated Weibull (EW) probability density function.
t: positive argument

sigma: scale parameter

nu: shape parameter

gamma: shape parameter

logd: log scale (TRUE or FALSE)
=#

function pdfEW(t, sigma, nu, gamma, logd::Bool=false)
    val = log(gamma) .+ (gamma - 1) * logcdf.(Weibull(nu, sigma), t) .+
          logpdf.(Weibull(nu, sigma), t)
    if logd
        return val
    else
        return exp.(val)
    end
end

#= Exponentiated Weibull (EW) cumulative distribution function.
t: positive argument

sigma: scale parameter

nu: shape parameter

gamma: shape parameter

logp: log scale (TRUE or FALSE)
=#

function cdfEW(t, sigma, nu, gamma, logp=false)
    val = gamma * logcdf.(Weibull(nu, sigma), t)
    if logp
        return val
    else
        return exp.(val)
    end
end


#= Exponentiated Weibull (EW) quantile function.
p: 	probability. A value in (0,1)

sigma: scale parameter

nu: shape parameter

gamma: shape parameter

logp: log scale (TRUE or FALSE)
=#

function qEW(p, sigma, nu, gamma)
    dist0 = Weibull(nu, sigma)
    prob = p .^ (1 / gamma)
    val = quantile(dist0, prob)
    return val
end

#= Exponentiated Weibull (EW) hazard function.
t: positive argument

sigma: scale parameter

nu: shape parameter

gamma: shape parameter

logh: log scale (TRUE or FALSE)
=#

function hEW(t, sigma, nu, gamma, logh::Bool=false)
    val = pdfEW.(t, sigma, nu, gamma, true) .-
          log.(1 .- cdfEW.(t, sigma, nu, gamma, false))
    if logh
        return val
    else
        return exp.(val)
    end
end

#= Exponentiated Weibull (EW) cumulative hazard function.
t: positive argument

sigma: scale parameter

nu: shape parameter

gamma: shape parameter
=#

function chEW(t, sigma, nu, gamma)
    val = -log.(1 .- cdfEW.(t, sigma, nu, gamma, false))
    return val
end


#= Exponentiated Weibull (EW) random number generation.
n: number of observations

sigma: scale parameter

nu: shape parameter

gamma: shape parameter
=#

function randEW(n::Int, sigma, nu, gamma)
    p = rand(n) .^ (1 / gamma)
    out = quantile.(Weibull(nu, sigma), p)
    return out
end