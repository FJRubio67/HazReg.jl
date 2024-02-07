"""
    simGH(...)


    Details...

Docs to rewrite....


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
beta  : regression parameters multiplying the hazard for the GH model
        or the regression parameters for AFT and PH models
alpha  : regression parameters multiplying the time scale for the GH model
        or the regression parameters for the AH model
des : Design matrix for the GH model (hazard scale)
      or design matrix for AFT and PH models
des_t : Design matrix for the GH model (time scale)
       or design matrix for the AH model
hstr  : hazard structure (AH, AFT, PH, GH)
dist: distribution passed as an already instanciated distribution from Distributions.jl :)
baseline  : baseline hazard distribution

Returns a vector containing the simulated times to event

References: 
* [Link to my reference so that people understand what it is](https://myref.com) 
"""
function simGH(; seed::Int64, n::Int64,
    des::Union{Matrix{Float64},Vector{Float64},Nothing}, 
    des_t::Union{Matrix{Float64},Vector{Float64},Nothing},
    alpha::Union{Vector{Float64},Float64,Nothing}, 
    beta::Union{Vector{Float64},Float64,Nothing},
    hstr::String, dist::Distributions.ContinuousDistribution)

    #= Uniform variates =#
    Random.seed!(seed)
    distu = Uniform(0, 1)
    u = rand(distu, n)
    d = dist # alias.

    # Seing this I am starting to think that here and in GMLE.jl these names should actually be structs, with a `loglikelyhood(...)` method for the computation done in GMLE, and a `predict(...)` method that does this thing. Would make the code much clearer IMHO. 

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
    times = quantile.(d,p0) ./ exp_xalpha

    return times
end