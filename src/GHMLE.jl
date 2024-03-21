#= 
=========================================================================
Optimisation function
========================================================================= 
=#

#= GHMLE function. 

Maximum likelihood estimation in General Hazards models using 
several parametric baseline hazards

init : initial point for optimisation step under the parameterisation 
(log(scale), log(shape1), log(shape2), alpha, beta) for scale-shape1-shape2 models or 
(mu, log(scale), alpha, beta) for log-location scale models.

times : times to event

status: vital status indicators (true or 1 = observed, false or 0 = censored)

hstr: hazard structure.  No covariates ("baseline"), AFT model ("AFT"), PH model ("PH"), 
      AH model ("AH"), GH model ("GH").

dist: baseline hazard distribution (LogNormal, LogLogistic, Weibull, Gamma, PGW, GenGamma, EW)

des: design matrix for hazard-level effects

des_t: design matrix for time-level effects (it is recommended not to use splines here)

method: "NM" (NelderMead), "N" (Newton), "LBFGS" (LBFGS), "CG" (ConjugateGradient), "GD" (GradientDescent).  

maxit: maximum number of iterations in "method"
=#

function GHMLE(; init::Vector{Float64}, times::Vector{Float64}, status::Vector{Bool},
    hstr::String, dist::String,
    des::Union{Matrix{Float64},Vector{Float64},Nothing}, 
    des_t::Union{Matrix{Float64},Vector{Float64},Nothing},
    method::String, maxit::Int64)

    #= -Log-likelihood =#
    function mloglik(par::Vector)
        #= 
        ****************************************************
        Baseline Hazards models 
        ****************************************************
        =#
        if hstr == "baseline"

            #= LogNormal baseline model =#
            if dist == "LogNormal"
                #= -Log-likelihood value =#
                mu = par[1]
                sigma = exp(par[2])
                val = -sum(hLogNormal(times[status], mu, sigma, true)) .+
                      sum(chLogNormal(times, mu, sigma))
            end

            #= LogLogistic baseline model =#
            if dist == "LogLogistic"
                #= -Log-likelihood value =#
                mu = par[1]
                sigma = exp(par[2])
                val = -sum(hLogLogistic(times[status], mu, sigma, true)) .+
                      sum(chLogLogistic(times, mu, sigma))
            end

            #= Weibull baseline model =#
            if dist == "Weibull"
                #= -Log-likelihood value =#
                shape = exp(par[1])
                scale = exp(par[2])
                val = -sum(hWeibull(times[status], shape, scale, true)) .+
                      sum(chWeibull(times, shape, scale))
            end

            #= Gamma baseline model =#
            if dist == "Gamma"
                #= -Log-likelihood value =#
                shape = exp(par[1])
                scale = exp(par[2])
                val = -sum(hGamma(times[status], shape, scale, true)) .+
                      sum(chGamma(times, shape, scale))
            end

            #= PGW baseline model =#
            if dist == "PGW"
                #= -Log-likelihood value =#
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                val = -sum(hPGW(times[status], sigma, nu, gamma, true)) .+
                      sum(chPGW(times, sigma, nu, gamma))
            end

            #= EW baseline model =#
            if dist == "EW"
                #= -Log-likelihood value =#
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                val = -sum(hEW(times[status], sigma, nu, gamma, true)) .+
                      sum(chEW(times, sigma, nu, gamma))
            end

            #= GenGamma baseline model =#
            if dist == "GenGamma"
                #= -Log-likelihood value =#
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                val = -sum(hGenGamma(times[status], sigma, nu, gamma, true)) .+
                      sum(chGenGamma(times, sigma, nu, gamma))
            end
        end

        #= 
        ****************************************************
        Proportional Hazards models 
        ****************************************************
        =#
        if hstr == "PH"

            #= LogNormal baseline model =#
            if dist == "LogNormal"
                #= -Log-likelihood value =#
                p = size(des, 2)
                mu = par[1]
                sigma = exp(par[2])
                beta = par[3:p+2]
                if p==1
                    beta = par[3]
                end

                x_beta = des * beta
                x_beta_obs = x_beta[status]
                exp_beta = exp.(x_beta)
                val = -sum(hLogNormal(times[status], mu, sigma, true) .+ x_beta_obs) .+
                      sum(chLogNormal(times, mu, sigma) .* exp_beta)
            end

            #= LogLogistic baseline model =#
            if dist == "LogLogistic"
                #= -Log-likelihood value =#
                p = size(des, 2)
                mu = par[1]
                sigma = exp(par[2])
                beta = par[3:p+2]
                if p==1
                    beta = par[3]
                end

                x_beta = des * beta
                x_beta_obs = x_beta[status]
                exp_beta = exp.(x_beta)
                val = -sum(hLogLogistic(times[status], mu, sigma, true) .+ x_beta_obs) .+
                      sum(chLogLogistic(times, mu, sigma) .* exp_beta)
            end

            #= Weibull baseline model =#
            if dist == "Weibull"
                #= -Log-likelihood value =#
                p = size(des, 2)
                shape = exp(par[1])
                scale = exp(par[2])
                beta = par[3:p+2]
                if p==1
                    beta = par[3]
                end

                x_beta = des * beta
                x_beta_obs = x_beta[status]
                exp_beta = exp.(x_beta)
                val = -sum(hWeibull(times[status], shape, scale, true) .+ x_beta_obs) .+
                      sum(chWeibull(times, shape, scale) .* exp_beta)
            end

            #= Gamma baseline model =#
            if dist == "Gamma"
                #= -Log-likelihood value =#
                p = size(des, 2)
                shape = exp(par[1])
                scale = exp(par[2])
                beta = par[3:p+2]
                if p==1
                    beta = par[3]
                end

                x_beta = des * beta
                x_beta_obs = x_beta[status]
                exp_beta = exp.(x_beta)
                val = -sum(hGamma(times[status], shape, scale, true) .+ x_beta_obs) .+
                      sum(chGamma(times, shape, scale) .* exp_beta)
            end

            #= PGW baseline model =#
            if dist == "PGW"
                #= -Log-likelihood value =#
                p = size(des, 2)
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                beta = par[4:p+3]
                if p==1
                    beta = par[4]
                end

                x_beta = des * beta
                x_beta_obs = x_beta[status]
                exp_beta = exp.(x_beta)
                val = -sum(hPGW(times[status], sigma, nu, gamma, true) .+ x_beta_obs) .+
                      sum(chPGW(times, sigma, nu, gamma) .* exp_beta)
            end

            #= EW baseline model =#
            if dist == "EW"
                #= -Log-likelihood value =#
                p = size(des, 2)
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                beta = par[4:p+3]
                if p==1
                    beta = par[4]
                end

                x_beta = des * beta
                x_beta_obs = x_beta[status]
                exp_beta = exp.(x_beta)
                val = -sum(hEW(times[status], sigma, nu, gamma, true) .+ x_beta_obs) .+
                      sum(chEW(times, sigma, nu, gamma) .* exp_beta)
            end

            #= GenGamma baseline model =#
            if dist == "GenGamma"
                #= -Log-likelihood value =#
                p = size(des, 2)
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                beta = par[4:p+3]
                if p==1
                    beta = par[4]
                end

                x_beta = des * beta
                x_beta_obs = x_beta[status]
                exp_beta = exp.(x_beta)
                val = -sum(hGenGamma(times[status], sigma, nu, gamma, true) .+ x_beta_obs) .+
                      sum(chGenGamma(times, sigma, nu, gamma) .* exp_beta)
            end
        end

        #= 
        ****************************************************
        Accelerated Failure Time models 
        ****************************************************
        =#
        if hstr == "AFT"

            #= LogNormal baseline model =#
            if dist == "LogNormal"
                #= -Log-likelihood value =#
                p = size(des, 2)
                mu = par[1]
                sigma = exp(par[2])
                beta = par[3:p+2]
                if p==1
                    beta = par[3]
                end

                x_beta = des * beta
                x_beta_obs = x_beta[status]
                exp_beta = exp.(x_beta)
                exp_beta_obs = exp.(x_beta_obs)
                val = -sum(hLogNormal(times[status] .* exp_beta_obs, mu, sigma, true) .+ x_beta_obs) .+
                      sum(chLogNormal(times .* exp_beta, mu, sigma))
            end

            #= LogLogistic baseline model =#
            if dist == "LogLogistic"
                #= -Log-likelihood value =#
                p = size(des, 2)
                mu = par[1]
                sigma = exp(par[2])
                beta = par[3:p+2]
                if p==1
                    beta = par[3]
                end

                x_beta = des * beta
                x_beta_obs = x_beta[status]
                exp_beta = exp.(x_beta)
                exp_beta_obs = exp.(x_beta_obs)
                val = -sum(hLogLogistic(times[status] .* exp_beta_obs, mu, sigma, true) .+ x_beta_obs) .+
                      sum(chLogLogistic(times .* exp_beta, mu, sigma))
            end

            #= Weibull baseline model =#
            if dist == "Weibull"
                #= -Log-likelihood value =#
                p = size(des, 2)
                shape = exp(par[1])
                scale = exp(par[2])
                beta = par[3:p+2]
                if p==1
                    beta = par[3]
                end

                x_beta = des * beta
                x_beta_obs = x_beta[status]
                exp_beta = exp.(x_beta)
                exp_beta_obs = exp.(x_beta_obs)
                val = -sum(hWeibull(times[status] .* exp_beta_obs, shape, scale, true) .+ x_beta_obs) .+
                      sum(chWeibull(times .* exp_beta, shape, scale))
            end

            #= Gamma baseline model =#
            if dist == "Gamma"
                #= -Log-likelihood value =#
                p = size(des, 2)
                shape = exp(par[1])
                scale = exp(par[2])
                beta = par[3:p+2]
                if p==1
                    beta = par[3]
                end

                x_beta = des * beta
                x_beta_obs = x_beta[status]
                exp_beta = exp.(x_beta)
                exp_beta_obs = exp.(x_beta_obs)
                val = -sum(hGamma(times[status] .* exp_beta_obs, shape, scale, true) .+ x_beta_obs) .+
                      sum(chGamma(times .* exp_beta, shape, scale))
            end

            #= PGW baseline model =#
            if dist == "PGW"
                #= -Log-likelihood value =#
                p = size(des, 2)
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                beta = par[4:p+3]
                if p==1
                    beta = par[4]
                end

                x_beta = des * beta
                x_beta_obs = x_beta[status]
                exp_beta = exp.(x_beta)
                exp_beta_obs = exp.(x_beta_obs)
                val = -sum(hPGW(times[status] .* exp_beta_obs, sigma, nu, gamma, true) .+ x_beta_obs) .+
                      sum(chPGW(times .* exp_beta, sigma, nu, gamma))
            end

            #= EW baseline model =#
            if dist == "EW"
                #= -Log-likelihood value =#
                p = size(des, 2)
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                beta = par[4:p+3]
                if p==1
                    beta = par[4]
                end

                x_beta = des * beta
                x_beta_obs = x_beta[status]
                exp_beta = exp.(x_beta)
                exp_beta_obs = exp.(x_beta_obs)
                val = -sum(hEW(times[status] .* exp_beta_obs, sigma, nu, gamma, true) .+ x_beta_obs) .+
                      sum(chEW(times .* exp_beta, sigma, nu, gamma))
            end

            #= GenGamma baseline model =#
            if dist == "GenGamma"
                #= -Log-likelihood value =#
                p = size(des, 2)
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                beta = par[4:p+3]
                if p==1
                    beta = par[4]
                end

                x_beta = des * beta
                x_beta_obs = x_beta[status]
                exp_beta = exp.(x_beta)
                exp_beta_obs = exp.(x_beta_obs)
                val = -sum(hGenGamma(times[status] .* exp_beta_obs, sigma, nu, gamma, true) .+ x_beta_obs) .+
                      sum(chGenGamma(times .* exp_beta, sigma, nu, gamma))
            end
        end

        #= 
        ****************************************************
        Accelerated Hazards models 
        ****************************************************
        =#
        if hstr == "AH"

            #= LogNormal baseline model =#
            if dist == "LogNormal"
                #= -Log-likelihood value =#
                q = size(des_t, 2)
                mu = par[1]
                sigma = exp(par[2])
                alpha = par[3:q+2]
                if q==1
                    alpha = par[3]
                end

                x_alpha = des_t * alpha
                x_alpha_obs = x_alpha[status]
                exp_alpha = exp.(x_alpha)
                exp_malpha = exp.(-x_alpha)
                exp_alpha_obs = exp.(x_alpha_obs)
                val = -sum(hLogNormal(times[status] .* exp_alpha_obs, mu, sigma, true)) .+
                      sum(chLogNormal(times .* exp_alpha, mu, sigma) .* exp_malpha)
            end

            #= LogLogistic baseline model =#
            if dist == "LogLogistic"
                #= -Log-likelihood value =#
                q = size(des_t, 2)
                mu = par[1]
                sigma = exp(par[2])
                alpha = par[3:q+2]
                if q==1
                    alpha = par[3]
                end

                x_alpha = des_t * alpha
                x_alpha_obs = x_alpha[status]
                exp_alpha = exp.(x_alpha)
                exp_malpha = exp.(-x_alpha)
                exp_alpha_obs = exp.(x_alpha_obs)
                val = -sum(hLogLogistic(times[status] .* exp_alpha_obs, mu, sigma, true)) .+
                      sum(chLogLogistic(times .* exp_alpha, mu, sigma) .* exp_malpha)
            end

            #= Weibull baseline model =#
            if dist == "Weibull"
                #= -Log-likelihood value =#
                q = size(des_t, 2)
                shape = exp(par[1])
                scale = exp(par[2])
                alpha = par[3:q+2]
                if q==1
                    alpha = par[3]
                end

                x_alpha = des_t * alpha
                x_alpha_obs = x_alpha[status]
                exp_alpha = exp.(x_alpha)
                exp_malpha = exp.(-x_alpha)
                exp_alpha_obs = exp.(x_alpha_obs)
                val = -sum(hWeibull(times[status] .* exp_alpha_obs, shape, scale, true)) .+
                      sum(chWeibull(times .* exp_alpha, shape, scale) .* exp_malpha)
            end

            #= Gamma baseline model =#
            if dist == "Gamma"
                #= -Log-likelihood value =#
                q = size(des_t, 2)
                shape = exp(par[1])
                scale = exp(par[2])
                alpha = par[3:q+2]
                if q==1
                    alpha = par[3]
                end

                x_alpha = des_t * alpha
                x_alpha_obs = x_alpha[status]
                exp_alpha = exp.(x_alpha)
                exp_malpha = exp.(-x_alpha)
                exp_alpha_obs = exp.(x_alpha_obs)
                val = -sum(hGamma(times[status] .* exp_alpha_obs, shape, scale, true)) .+
                      sum(chGamma(times .* exp_alpha, shape, scale) .* exp_malpha)
            end

            #= PGW baseline model =#
            if dist == "PGW"
                #= -Log-likelihood value =#
                q = size(des_t, 2)
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                alpha = par[4:q+3]
                if q==1
                    alpha = par[4]
                end

                x_alpha = des_t * alpha
                x_alpha_obs = x_alpha[status]
                exp_alpha = exp.(x_alpha)
                exp_malpha = exp.(-x_alpha)
                exp_alpha_obs = exp.(x_alpha_obs)
                val = -sum(hPGW(times[status] .* exp_alpha_obs, sigma, nu, gamma, true)) .+
                      sum(chPGW(times .* exp_alpha, sigma, nu, gamma) .* exp_malpha)
            end

            #= EW baseline model =#
            if dist == "EW"
                #= -Log-likelihood value =#
                q = size(des_t, 2)
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                alpha = par[4:q+3]
                if q==1
                    alpha = par[4]
                end

                x_alpha = des_t * alpha
                x_alpha_obs = x_alpha[status]
                exp_alpha = exp.(x_alpha)
                exp_malpha = exp.(-x_alpha)
                exp_alpha_obs = exp.(x_alpha_obs)
                val = -sum(hEW(times[status] .* exp_alpha_obs, sigma, nu, gamma, true)) .+
                      sum(chEW(times .* exp_alpha, sigma, nu, gamma) .* exp_malpha)
            end

            #= GenGamma baseline model =#
            if dist == "GenGamma"
                #= -Log-likelihood value =#
                q = size(des_t, 2)
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                alpha = par[4:q+3]
                if q==1
                    alpha = par[4]
                end

                x_alpha = des_t * alpha
                x_alpha_obs = x_alpha[status]
                exp_alpha = exp.(x_alpha)
                exp_malpha = exp.(-x_alpha)
                exp_alpha_obs = exp.(x_alpha_obs)
                val = -sum(hGenGamma(times[status] .* exp_alpha_obs, sigma, nu, gamma, true)) .+
                      sum(chGenGamma(times .* exp_alpha, sigma, nu, gamma) .* exp_malpha)
            end
        end

        #= 
        ****************************************************
        General Hazards models 
        ****************************************************
        =#
        if hstr == "GH"

            #= LogNormal baseline model =#
            if dist == "LogNormal"
                #= -Log-likelihood value =#
                p = size(des, 2)
                q = size(des_t, 2)
                mu = par[1]
                sigma = exp(par[2])
                alpha = par[3:q+2]
                beta = par[q+3:p+q+2]
                if q==1
                    alpha = par[3]
                end
                if p==1
                    beta = par[q+3]
                end

                x_alpha = des_t * alpha
                x_beta = des * beta
                x_alpha_obs = x_alpha[status]
                x_beta_obs = x_beta[status]
                exp_alpha = exp.(x_alpha)
                exp_alpha_obs = exp.(x_alpha_obs)
                exp_beta = exp.(x_beta)
                exp_beta_obs = exp.(x_beta_obs)
                exp_dif = exp.(x_beta .- x_alpha)
                val = -sum(hLogNormal(times[status] .* exp_alpha_obs, mu, sigma, true) .+ x_beta_obs) .+
                      sum(chLogNormal(times .* exp_alpha, mu, sigma) .* exp_dif)
            end

            #= LogLogistic baseline model =#
            if dist == "LogLogistic"
                #= -Log-likelihood value =#
                p = size(des, 2)
                q = size(des_t, 2)
                mu = par[1]
                sigma = exp(par[2])
                alpha = par[3:q+2]
                beta = par[q+3:p+q+2]
                if q==1
                    alpha = par[3]
                end
                if p==1
                    beta = par[q+3]
                end

                x_alpha = des_t * alpha
                x_beta = des * beta
                x_alpha_obs = x_alpha[status]
                x_beta_obs = x_beta[status]
                exp_alpha = exp.(x_alpha)
                exp_alpha_obs = exp.(x_alpha_obs)
                exp_beta = exp.(x_beta)
                exp_beta_obs = exp.(x_beta_obs)
                exp_dif = exp.(x_beta .- x_alpha)
                val = -sum(hLogLogistic(times[status] .* exp_alpha_obs, mu, sigma, true) .+ x_beta_obs) .+
                      sum(chLogLogistic(times .* exp_alpha, mu, sigma) .* exp_dif)
            end

            #= Weibull baseline model =#
            if dist == "Weibull"
                #= -Log-likelihood value =#
                p = size(des, 2)
                q = size(des_t, 2)
                shape = exp(par[1])
                scale = exp(par[2])
                alpha = par[3:q+2]
                beta = par[q+3:p+q+2]
                if q==1
                    alpha = par[3]
                end
                if p==1
                    beta = par[q+3]
                end

                x_alpha = des_t * alpha
                x_beta = des * beta
                x_alpha_obs = x_alpha[status]
                x_beta_obs = x_beta[status]
                exp_alpha = exp.(x_alpha)
                exp_alpha_obs = exp.(x_alpha_obs)
                exp_beta = exp.(x_beta)
                exp_beta_obs = exp.(x_beta_obs)
                exp_dif = exp.(x_beta .- x_alpha)
                val = -sum(hWeibull(times[status] .* exp_alpha_obs, shape, scale, true) .+ x_beta_obs) .+
                      sum(chWeibull(times .* exp_alpha, shape, scale) .* exp_dif)
            end

            #= Gamma baseline model =#
            if dist == "Gamma"
                #= -Log-likelihood value =#
                p = size(des, 2)
                q = size(des_t, 2)
                shape = exp(par[1])
                scale = exp(par[2])
                alpha = par[3:q+2]
                beta = par[q+3:p+q+2]
                if q==1
                    alpha = par[3]
                end
                if p==1
                    beta = par[q+3]
                end
                
                x_alpha = des_t * alpha
                x_beta = des * beta
                x_alpha_obs = x_alpha[status]
                x_beta_obs = x_beta[status]
                exp_alpha = exp.(x_alpha)
                exp_alpha_obs = exp.(x_alpha_obs)
                exp_beta = exp.(x_beta)
                exp_beta_obs = exp.(x_beta_obs)
                exp_dif = exp.(x_beta .- x_alpha)
                val = -sum(hGamma(times[status] .* exp_alpha_obs, shape, scale, true) .+ x_beta_obs) .+
                      sum(chGamma(times .* exp_alpha, shape, scale) .* exp_dif)
            end

            #= PGW baseline model =#
            if dist == "PGW"
                #= -Log-likelihood value =#
                p = size(des, 2)
                q = size(des_t, 2)
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                alpha = par[4:q+3]
                beta = par[q+4:p+q+3]
                if q==1
                    alpha = par[4]
                end
                if p==1
                    beta = par[q+4]
                end
                
                x_alpha = des_t * alpha
                x_beta = des * beta
                x_alpha_obs = x_alpha[status]
                x_beta_obs = x_beta[status]
                exp_alpha = exp.(x_alpha)
                exp_alpha_obs = exp.(x_alpha_obs)
                exp_beta = exp.(x_beta)
                exp_beta_obs = exp.(x_beta_obs)
                exp_dif = exp.(x_beta .- x_alpha)
                val = -sum(hPGW(times[status] .* exp_alpha_obs, sigma, nu, gamma, true) .+ x_beta_obs) .+
                      sum(chPGW(times .* exp_alpha, sigma, nu, gamma) .* exp_dif)
            end

            #= EW baseline model =#
            if dist == "EW"
                #= -Log-likelihood value =#
                p = size(des, 2)
                q = size(des_t, 2)
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                alpha = par[4:q+3]
                beta = par[q+4:p+q+3]
                if q==1
                    alpha = par[4]
                end
                if p==1
                    beta = par[q+4]
                end

                x_alpha = des_t * alpha
                x_beta = des * beta
                x_alpha_obs = x_alpha[status]
                x_beta_obs = x_beta[status]
                exp_alpha = exp.(x_alpha)
                exp_alpha_obs = exp.(x_alpha_obs)
                exp_beta = exp.(x_beta)
                exp_beta_obs = exp.(x_beta_obs)
                exp_dif = exp.(x_beta .- x_alpha)
                val = -sum(hEW(times[status] .* exp_alpha_obs, sigma, nu, gamma, true) .+ x_beta_obs) .+
                      sum(chEW(times .* exp_alpha, sigma, nu, gamma) .* exp_dif)
            end

            #= GenGamma baseline model =#
            if dist == "GenGamma"
                #= -Log-likelihood value =#
                p = size(des, 2)
                q = size(des_t, 2)
                sigma = exp(par[1])
                nu = exp(par[2])
                gamma = exp(par[3])
                alpha = par[4:q+3]
                beta = par[q+4:p+q+3]
                if q==1
                    alpha = par[4]
                end
                if p==1
                    beta = par[q+4]
                end

                x_alpha = des_t * alpha
                x_beta = des * beta
                x_alpha_obs = x_alpha[status]
                x_beta_obs = x_beta[status]
                exp_alpha = exp.(x_alpha)
                exp_alpha_obs = exp.(x_alpha_obs)
                exp_beta = exp.(x_beta)
                exp_beta_obs = exp.(x_beta_obs)
                exp_dif = exp.(x_beta .- x_alpha)
                val = -sum(hGenGamma(times[status] .* exp_alpha_obs, sigma, nu, gamma, true) .+ x_beta_obs) .+
                      sum(chGenGamma(times .* exp_alpha, sigma, nu, gamma) .* exp_dif)
            end
        end

        #= return -Log-likelihood value =#
        return val
    end

    #= Optimisation step =#
    if method == "NM"
        optimiser = optimize(mloglik, init, method=NelderMead(), iterations=maxit)
    end
    if method == "N"
        optimiser = optimize(mloglik, init, method=Newton(), iterations=maxit)
    end
    if method == "LBFGS"
        optimiser = optimize(mloglik, init, method=LBFGS(), iterations=maxit)
    end
    if method == "CG"
        optimiser = optimize(mloglik, init, method=ConjugateGradient(), iterations=maxit)
    end
    if method == "GD"
        optimiser = optimize(mloglik, init, method=GradientDescent(), iterations=maxit)
    end

    #= Returns the negative log-likelihood and the optimisation result =#
    return optimiser, mloglik
end

