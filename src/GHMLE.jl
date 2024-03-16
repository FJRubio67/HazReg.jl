"""
    GHMLE(...)


    Maximum likelihood estimation in General Hazards models using 
    several parametric baseline hazards

Docs to rewrite....


init : initial point for optimisation step under the parameterisation 
(log(scale), log(shape1), log(shape2), alpha, beta) for scale-shape1-shape2 models or 
(mu, log(scale), alpha, beta) for log-location scale models.

times : times to event

status: vital status indicators (true or 1 = observed, false or 0 = censored)

hstr: hazard structure.  No covariates ("baseline"), AFT model ("AFT"), PH model ("PH"), 
        AH model ("AH"), GH model ("GH").

dist: baseline hazard distribution, on the form of an already instanciated Distributions.jl's ContinuousUnivariateDistribution. Tested with LogNormal(), LogLogistic(), Weibull(), Gamma(), PGW(), EW(), GenGamma()

des: design matrix for hazard-level effects

des_t: design matrix for time-level effects (it is recommended not to use splines here)

method: one of NelderMead(), Newton(), LBFGS(), ConjugateGradient() or GradientDescent() or any other method taken by Optim.optimize()... 

maxit: maximum number of iterations of the optimization routine. 

References: 
* [Link to my reference so that people understand what it is](https://myref.com) 
"""
function GHMLE(; 
    init::Vector{Float64}, 
    times::Vector{Float64}, 
    status::Vector{Bool},
    hstr::String, 
    dist::Distributions.ContinuousUnivariateDistribution,
    des::Union{Matrix{Float64},Vector{Float64},Nothing}, 
    des_t::Union{Matrix{Float64},Vector{Float64},Nothing},
    method::Optim.AbstractOptimizer, # Use with method =  directly ! 
    maxit::Int64)
    
    # Quick side note : over-typing stuff in the function declaration only constraint your usage of the function and dos NOT add performance at all. This is only type restriction, the compiler will always compile on the acrtual types of the passed arguments regardless of what you input here. 
    # So you better NOT used these... 

    #= -Log-likelihood =#

    # fix nothings be replacing by empty matrices. 
    des = isnothing(des) ? Float64[;;] : des
    des_t = isnothing(des_t) ? Float64[;;] : des_t

    # get number of parameters: 
    npd = length(params(d))
    p = size(des,2)
    q = size(des_t, 2)
    @assert length(init) == npd + p + q # check init vector lenght: 

    function mloglik(par::Vector)
        # split parameters: 
        d  = typeof(dist)(par[1:npd]...) # reconstruct the distribution with these new parameters. 

        if (q == 0) && (p == 0)
            # then only the baseline would be returned since there is no parameters. 
            return - sum(loghaz.(d,time[status])) + sum(cumhaz.(d,times))
        end

        if p > 0 # then you might want some regression parameters beta. 
            beta = par[(npd+q+1):(npd+q+p)]
            x_beta = des * beta
            x_beta_obs = x_beta[status]
            exp_beta = exp.(x_beta)
        end

        if q > 0 # The you might ALSO want some parameters alpha. 
            alpha = par[(npd+1):(npd+q)]
            x_alpha = des_t * alpha
            x_alpha_obs = x_alpha[status]
            exp_alpha = exp.(x_alpha)
            exp_malpha = exp.(-x_alpha)
            exp_alpha_obs = exp.(x_alpha_obs)
        end

        # Proportional Hazards models 
        if hstr == "PH" 
            return -sum(loghaz.(d,times[status])) + sum(x_beta_obs) + sum(cumhaz.(d,times) .* exp_beta)
        end

        # Accelerated Failure Time models 
        if hstr == "AFT"
            exp_beta_obs = exp.(x_beta_obs)
            return -sum(loghaz.(d,times[status] .* exp_beta_obs)) + sum(x_beta_obs) + sum(cumhaz.(d,times .* exp_beta))
        end

        # Accelerated Hazards models 
        if hstr == "AH"
            return -sum(loghaz.(d,times[status] .* exp_alpha_obs)) + sum(cumhaz(d,times .* exp_alpha) .* exp_malpha)
        end

        # General Hazards models 
        if hstr == "GH"
            exp_dif = exp.(x_beta .- x_alpha)
            return -sum(loghaz.(d,times[status] .* exp_alpha_obs) .+ x_beta_obs) + sum(cumhaz(d, times .* exp_alpha) .* exp_dif)
        end

        # if you end up here then you have provided wrong arguments. 
        @error("Wrong arguments.")
    end

    optimiser = optimize(mloglik, init, method=NelderMead(), iterations=maxit)
    return optimiser, mloglik
end


