#= 
----------------------------------------------------------------------------------------
Function to calculate normal confidence intervals.
It uses the reparameterised log-likelihood function, where all positive parameters are
mapped to the real line using a log-link.
----------------------------------------------------------------------------------------
=#

#= 
FUN   : minus log-likelihood function to be used to calculate the confidence intervals
MLE   : maximum likelihood estimator of the parameters 
level : confidence level
Returns a list containing the upper and lower conf.int limits, the transformed MLE, and std errors
=#


function ConfInt(; FUN::Function, MLE::Vector{Float64}, level::Float64)
    dist = Normal(0,1)
    sd_int = abs(quantile(dist, 0.5*(1-level)))

    HESS = ForwardDiff.hessian(FUN, MLE)
    Fisher_Info = inv(HESS)
    Sigma = sqrt.(diag(Fisher_Info))
    U = MLE .+ sd_int .* Sigma
    L = MLE .- sd_int .* Sigma
    CI = hcat(L,U)
    
    return CI
end