#=
Power Generalised Weibull Distribution: Hazard, cumulative hazard,
probability density function, random number generation, and survival function
=#

#= Power Generalised Weibull (PGW) hazard function.
t: positive argument

sigma: scale parameter

nu: shape parameter

gamma: shape parameter

logh: log scale (TRUE or FALSE)
=#
function hPGW(t, sigma, nu, gamma, logh::Bool=false)
    val = log.(nu) .- log.(gamma) .- nu * log.(sigma) .+ (nu - 1) * log.(t) .+
          (1 / gamma - 1) * log.(1 .+ (t / sigma) .^ nu)
    if logh
        return val
    else
        return exp.(val)
    end
end

#= Power Generalised Weibull (PGW) cumulative hazard function.
t: positive argument

sigma: scale parameter

nu: shape parameter

gamma: shape parameter
=#
function chPGW(t, sigma, nu, gamma)
    val = -1 .+ (1 .+ (t / sigma) .^ nu) .^ (1 / gamma)
    return val
end

#= Power Generalised Weibull (PGW) probability density function.
t: positive argument

sigma: scale parameter

nu: shape parameter

gamma: shape parameter

logd: log scale (TRUE or FALSE)
=#
function pdfPGW(t, sigma, nu, gamma, logd::Bool=false)
    val = log.(nu) .- log.(gamma) .- nu * log.(sigma) .+ (nu .- 1) * log.(t) .+
          (1 / gamma - 1) * log.(1 .+ (t / sigma) .^ nu) .+
          (1 .- (1 .+ (t / sigma) .^ nu) .^ (1 / gamma))
    if logd
        return val
    else
        return exp.(val)
    end
end


#= Power Generalised Weibull (PGW) survival function.
t: positive argument

sigma: scale parameter

nu: shape parameter

gamma: shape parameter

logp: log scale (TRUE or FALSE)
=#

function ccdfPGW(t, sigma, nu, gamma, logp=false)
    val = 1 .- (1 .+ (t / sigma) .^ nu) .^ (1 / gamma)
    if logp
        return val
    else
        return exp.(val)
    end
end

#= Power Generalised Weibull (PGW) quantile function.
p: 	probability. A value in (0,1)

sigma: scale parameter

nu: shape parameter

gamma: shape parameter
=#
function qPGW(p, sigma, nu, gamma)
    val = sigma * ((1 .- log.(1 .- p)) .^ gamma .- 1) .^ (1 / nu)
    return val
end



#= Power Generalised Weibull (PGW) random number generation.
n: number of observations

sigma: scale parameter

nu: shape parameter

gamma: shape parameter
=#

function randPGW(n::Int, sigma, nu, gamma)
    p = rand(n)
    out = sigma * ((1 .- log.(1 .- p)) .^ gamma .- 1) .^ (1 / nu)
    return out
end