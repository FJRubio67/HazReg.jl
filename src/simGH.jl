#= 
--------------------------------------------------------------------------------------------------
simGH function: Function to simulate times to event from a model with AH, AFT, PH, GH structures
for different parametric baseline hazards.
Distributions: LogNormal, LogLogistic, GenGamma, GGamma, Weibull, PGW, EW.
See: https://github.com/FJRubio67/HazReg
--------------------------------------------------------------------------------------------------
=#

#=
seed  : seed for simulation
n : sample size (number of individuals)
theta  :  parameters of the baseline hazard
beta  : regression parameters multiplying the hazard for the GH model
        or the regression parameters for AFT and PH models
alpha  : regression parameters multiplying the time scale for the GH model
        or the regression parameters for the AH model
des : Design matrix for the GH model (hazard scale)
      or design matrix for AFT and PH models
des_t : Design matrix for the GH model (time scale)
       or design matrix for the AH model
hstr  : hazard structure (AH, AFT, PH, GH)
baseline  : baseline hazard distribution

Returns a vector containing the simulated times to event
=#

function simGH(; seed::Int64, n::Int64,
    des::Union{Matrix{Float64},Vector{Float64},Nothing}, 
    des_t::Union{Matrix{Float64},Vector{Float64},Nothing},
    theta::Vector{Float64}, 
    alpha::Union{Vector{Float64},Float64,Nothing}, 
    beta::Union{Vector{Float64},Float64,Nothing},
    hstr::String, dist::String)

    #= Uniform variates =#
    Random.seed!(seed)
    distu = Uniform(0, 1)
    u = rand(distu, n)

    #= quantile function =#
    function quantf(prob)

        #= LogNormal baseline model =#
        if dist == "LogNormal"
            #= quantile value =#
            mu = theta[1]
            sigma = theta[2]
            distLN = LogNormal(mu, sigma)
            val = quantile(distLN, prob)
        end

        #= LogLogistic baseline model =#
        if dist == "LogLogistic"
            #= quantile value =#
            mu = theta[1]
            sigma = theta[2]
            distLL = LogLogistic(mu, sigma)
            val = quantile(distLL, prob)
        end

        #= Weibull baseline model =#
        if dist == "Weibull"
            #= quantile value =#
            shape = theta[1]
            scale = theta[2]
            distW = Weibull(shape, scale)
            val = quantile(distW, prob)
        end

        #= Gamma baseline model =#
        if dist == "Gamma"
            #= quantile value =#
            shape = theta[1]
            scale = theta[2]
            distG = Gamma(shape, scale)
            val = quantile(distG, prob)
        end

        #= PGW baseline model =#
        if dist == "PGW"
            #= quantile value =#
            sigma = theta[1]
            nu = theta[2]
            gamma = theta[3]
            val = qPGW(prob, sigma, nu, gamma)
        end

        #= GenGamma baseline model =#
        if dist == "EW"
            #= quantile value =#
            sigma = theta[1]
            nu = theta[2]
            gamma = theta[3]
            val = qEW(prob, sigma, nu, gamma)
        end

        #= GenGamma baseline model =#
        if dist == "GenGamma"
            #= quantile value =#
            sigma = theta[1]
            nu = theta[2]
            gamma = theta[3]
            val = qGenGamma(prob, sigma, nu, gamma)
        end

        return val
    end

    #= GH model =#
    if hstr == "GH"
        # Linear predictors
        exp_xalpha = exp.(des_t * alpha)
        exp_dif = exp.(des_t * alpha .- des * beta)
    end

    #= PH model =#
    if hstr == "PH"
        # Linear predictors
        exp_xalpha = 1.0
        exp_dif = exp.(-des * beta)
    end

    #= AFT model =#
    if hstr == "AFT"
        # Linear predictors
        exp_xalpha = exp.(des * beta)
        exp_dif = 1.0
    end

    #= AH model =#
    if hstr == "AH"
        # Linear predictors
        exp_xalpha = exp.(des_t * alpha)
        exp_dif = exp.(des_t * alpha)
    end

    # Simulating the times to event
    p0 = 1 .- exp.(log.(1 .- u) .* exp_dif)
    times = quantf.(p0) ./ exp_xalpha

    return times
end