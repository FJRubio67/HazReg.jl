#=
Generalised Gamma (GG) Distribution: Hazard, cumulative hazard,
probability density function, random number generation, and survival function
=#


#= Generalised Gamma (GG) probability density function.
t: positive argument

sigma: scale parameter

nu: shape parameter

gamma: shape parameter

logd: log scale (TRUE or FALSE)
=#

function pdfGenGamma(t, sigma, nu, gamma, logd::Bool=false)
    val = log(gamma) .- nu * log(sigma) .- loggamma(nu / gamma) .+
          (nu - 1) * log.(t) .- (t / sigma) .^ gamma
    if logd
        return val
    else
        return exp.(val)
    end
end

#= Generalised Gamma (GG) survival function.
t: positive argument

sigma: scale parameter

nu: shape parameter

gamma: shape parameter

logp: log scale (TRUE or FALSE)
=#

function ccdfGenGamma(t, sigma, nu, gamma, logp::Bool=false)
    tp = t .^ gamma
    val = logccdf.(Gamma(nu / gamma, sigma^gamma), tp)
    if logp
        return val
    else
        return exp.(val)
    end
end

#= Generalised Gamma (GG) quantile function.
p: 	probability. A value in (0,1)

sigma: scale parameter

nu: shape parameter

gamma: shape parameter

logp: log scale (TRUE or FALSE)
=#

function qGenGamma(p, sigma, nu, gamma)
    dist0 = Gamma(nu / gamma, sigma^gamma)
    val = quantile(dist0, p) .^ (1 / gamma)
    return val
end

#= Generalised Gamma (GG) hazard function.
t: positive argument

sigma: scale parameter

nu: shape parameter

gamma: shape parameter

logh: log scale (TRUE or FALSE)
=#

function hGenGamma(t, sigma, nu, gamma, logh::Bool=false)
    lpdf0 = pdfGenGamma.(t, sigma, nu, gamma, true)
    ls0 = ccdfGenGamma.(t, sigma, nu, gamma, true)
    val = lpdf0 .- ls0
    if logh
        return val
    else
        return exp.(val)
    end
end

#= Generalised Gamma (GG) cumulative hazard function.
t: positive argument

sigma: scale parameter

nu: shape parameter

gamma: shape parameter
=#

function chGenGamma(t, sigma, nu, gamma)
    H0 = -ccdfGenGamma.(t, sigma, nu, gamma, true)
    return H0
end


#= Generalised Gamma (GG) random number generation.
n: number of observations

sigma: scale parameter

nu: shape parameter

gamma: shape parameter
=#

function randGenGamma(n::Int, sigma, nu, gamma)
    p = rand(n)
    out = quantile.(Gamma(nu / gamma, sigma^gamma), p) .^ (1 / gamma)
    return out
end