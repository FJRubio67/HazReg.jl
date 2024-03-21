# HazReg.jl: Parametric Hazard-based regression models for survival data



# Models

The `HazReg.jl` Julia package implements the following parametric hazard-based
regression models for (overall) survival data.

-   General Hazard (GH) model [chen:2001](@cite) [rubio:2019](@cite).

-   Accelerated Failure Time (AFT) model [kalbfleisch:2011](@cite).

-   Proportional Hazards (PH) model [cox:1972](@cite).

-   Accelerated Hazards (AH) model [chen:2000](@cite).

These models are fitted using the Julia package `Optim` (methods included: "NM" (NelderMead), "N" (Newton), "LBFGS" (LBFGS), "CG" (ConjugateGradient), "GD" (GradientDescent)). Thus, the user needs to specify the initial points and to check the convergence of the optimisation step, as usual.

A description of these hazard models is presented below as well as the available baseline hazards.

## General Hazard model

The GH model is formulated in terms of the hazard structure 

$$h(t; \alpha, \beta, \theta, {\bf x}) = h_0\left(t  \exp\{\tilde{\bf x}^{\top}\alpha\}; \theta\right) \exp\{{\bf x}^{\top}\beta\}.$$ 

where ${\bf x}\in{\mathbb R}^p$ are the covariates that affect the hazard level; $\tilde{\bf x} \in {\mathbb R}^q$ are the covariates the affect the time level (typically $\tilde{\bf x} \subset {\bf x}$); $\alpha \in {\mathbb R}^q$ and $\beta \in {\mathbb R}^p$ are the regression coefficients; and $\theta \in \Theta$ is the vector of parameters of the baseline hazard $h_0(\cdot)$.

This hazard structure leads to an identifiable model as long as the baseline hazard is not a hazard associated to a member of the Weibull family of distributions [chen:2001](@cite).

## Accelerated Failure Time (AFT) model

The AFT model is formulated in terms of the hazard structure 

$$h(t; \beta, \theta, {\bf x}) = h_0\left(t  \exp\{{\bf x}^{\top}\beta\}; \theta\right) \exp\{{\bf x}^{\top}\beta\}.$$

where ${\bf x}\in{\mathbb R}^p$ are the available covariates; $\beta \in {\mathbb R}^p$ are the regression coefficients; and
$\theta \in \Theta$ is the vector of parameters of the baseline hazard $h_0(\cdot)$.

## Proportional Hazards (PH) model

The PH model is formulated in terms of the hazard structure 

$$h(t; \beta, \theta, {\bf x}) = h_0\left(t ; \theta\right) \exp\{{\bf x}^{\top}\beta\}.$$

where ${\bf x}\in{\mathbb R}^p$ are the available covariates; $\beta \in {\mathbb R}^p$ are the regression coefficients; and
$\theta \in \Theta$ is the vector of parameters of the baseline hazard $h_0(\cdot)$.

## Accelerated Hazards (AH) model

The AH model is formulated in terms of the hazard structure 

$$h(t; \alpha, \theta, \tilde{\bf x}) = h_0\left(t \exp\{\tilde{\bf x}^{\top}\alpha\}; \theta\right) .$$

where $\tilde{\bf x}\in{\mathbb R}^q$ are the available covariates; $\alpha \in {\mathbb R}^q$ are the regression coefficients; and
$\theta \in \Theta$ is the vector of parameters of the baseline hazard $h_0(\cdot)$.

# Available baseline hazards

The current version of the `HazReg.jl` Julia package implements the following
parametric baseline hazards for the models discussed in the previous
section.

-   [Power Generalised Weibull](http://rpubs.com/FJRubio/PGW) (`PGW`)
    distribution.

-   [Exponentiated Weibull](http://rpubs.com/FJRubio/EWD) (`EW`)
    distribution.

-   [Generalised Gamma](http://rpubs.com/FJRubio/GG) (`GenGamma`) distribution.

-   [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution) (`Gamma`)
    distribution.

-   [Lognormal](https://en.wikipedia.org/wiki/Log-normal_distribution)
    (`LogNormal`) distribution.

-   [Log-logistic](https://en.wikipedia.org/wiki/Log-logistic_distribution)
    (`LogLogistic`) distribution.

-   [Weibull](https://en.wikipedia.org/wiki/Weibull_distribution) (`Weibull`)
    distribution. (only for AFT, PH, and AH models)

All positive parameters are transformed into the real line using a `log` link (reparameterisation).

# Illustrative example: Julia code

In this example, we analyse the `LeukSurv` data set from the R package
`spBayesSurv`. This data set contains information about the survival of
acute myeloid leukemia in 1,043 patients.

For the GH model, we consider the hazard level covariates (${\bf x}$)
age (standardised), sex, wbc (white blood cell count at diagnosis,
standardised), and tpi (the Townsend score, standardised); and the time
level covariates (${\bf x}$) age (standardised), wbc (white blood cell
count at diagnosis, standardised), and tpi (the Townsend score,
standardised). For the PH, AFT, and AH models, we consider the
covariates age (standardised), sex, wbc (white blood cell count at
diagnosis, standardised), and tpi (the Townsend score, standardised).

For illustration, we fit the 4 models with both (3-parameter) PGW and
(2-parameter) LL baseline hazard. In addition, we fit the GH model with
`GenGamma`, `EW`, `LogNormal`, `LogLogistic`, and `Gamma` baseline hazards. We compare these models in terms of
AIC (BIC can be used as well). We summarise the best selected model with
the available tools in this package.

See also: 

- [HazReg.jl Julia Package](https://github.com/FJRubio67/HazReg.jl) 

- [HazReg](https://github.com/FJRubio67/GHSurv), for an R implementation.


## Data preparation

```@example 1

using Distributions
using Random
using StatsBase
using Optim
using LinearAlgebra
using SpecialFunctions
using ForwardDiff


using HazReg
using Plots
using PrettyTables
using DataFrames
using NamedArrays
using CSV
using Survival

#= Data =#
df = CSV.File(joinpath(@__DIR__,"..","src","assets","LeukSurv.csv"));

#= Design matrix for hazard level effects =#
des = hcat( standardise(df.age), df.sex, standardise(df.wbc), standardise(df.tpi) );

#= Design matrix for time level effects =#
des_t = hcat( standardise(df.age), standardise(df.wbc), standardise(df.tpi) );

#= Vital status =#
status = collect(Bool,(df.cens));

#= Survival times =#
times = df.time/365.25 ;
```

## Model fit and MLEs

```@example 1
# PGWGH
OPTPGWGH = GHMLE(init = fill(0.0, 3 + size(des_t)[2] + size(des)[2]), times = times,
            status = status, hstr = "GH", dist = "PGW", 
            des = des, des_t = des_t, method = "NM", maxit = 1000)

# PGWAFT
OPTPGWAFT = GHMLE(init = fill(0.0, 3 + size(des)[2]), times = times, 
                  status = status, hstr = "AFT", dist = "PGW", 
                  des = des, des_t = nothing, method = "NM", maxit = 1000)

# PGWPH
OPTPGWPH = GHMLE(init = fill(0.0, 3 +size(des)[2]), times = times, 
                 status = status, hstr = "PH", dist = "PGW", 
                 des = des, des_t = nothing, method = "NM", maxit = 1000)

# PGWAH
OPTPGWAH = GHMLE(init = fill(0.0, 3 + size(des_t)[2] ), times = times, 
                 status = status, hstr = "AH", dist = "PGW", 
                 des_t = des_t, des = nothing, method = "NM", maxit = 1000)


# LLGH
OPTLLGH = GHMLE(init = fill(0.0, 2 + size(des_t)[2] + size(des)[2]), times = times,
                status = status, hstr = "GH", dist = "LogLogistic", 
                des = des, des_t = des_t, method = "NM", maxit = 1000)

# LLAFT
OPTLLAFT = GHMLE(init = fill(0.0, 2 + size(des)[2]), times = times, 
                  status = status, hstr = "AFT", dist = "LogLogistic", 
                  des = des, des_t = nothing, method = "NM", maxit = 1000)

# LLPH
OPTLLPH = GHMLE(init = fill(0.0, 2 + size(des)[2]), times = times, 
                status = status, hstr = "PH", dist = "LogLogistic", 
                des = des, des_t = nothing, method = "NM", maxit = 1000)

# LLAH
OPTLLAH = GHMLE(init = fill(0.0, 2 + size(des_t)[2]), times = times, 
                status = status, hstr = "AH", dist = "LogLogistic", 
                des = nothing, des_t = des_t, method = "NM", maxit = 1000)


# EWGH
OPTEWGH = GHMLE(init = fill(0.0, 3 + size(des_t)[2] + size(des)[2]), times = times, 
                status = status, hstr = "GH", dist = "EW", 
                des = des, des_t = des_t, method = "NM", maxit = 1000)

# GGGH
OPTGGGH = GHMLE(init = fill(0.0, 3 + size(des_t)[2] + size(des)[2]), times = times,
                status = status, hstr = "GH", dist = "GenGamma",
                des = des, des_t = des_t, method = "NM", maxit = 1000)

# LNGH
OPTLNGH = GHMLE(init = fill(0.0, 2 + size(des_t)[2] + size(des)[2]), times = times, 
                status = status, hstr = "GH", dist = "LogNormal", 
                des = des, des_t = des_t, method = "NM", maxit = 1000)

# GGH
OPTGGH = GHMLE(init = fill(0.0, 2 + size(des_t)[2] + size(des)[2]), times = times, 
               status = status, hstr = "GH", dist = "Gamma", 
               des = des, des_t = des_t, method = "N", maxit = 1000)

# MLEs in the original parameterisations
MLEPGWGH = [exp(OPTPGWGH[1].minimizer[j]) for j in 1:3] 
append!(MLEPGWGH, OPTPGWGH[1].minimizer[4:end]) 

MLEEWGH = [exp(OPTEWGH[1].minimizer[j]) for j in 1:3] 
append!(MLEEWGH, OPTEWGH[1].minimizer[4:end]) 

MLEEWGH = [exp(OPTEWGH[1].minimizer[j]) for j in 1:3] 
append!(MLEEWGH, OPTEWGH[1].minimizer[4:end]) 

MLEGGGH = [exp(OPTGGGH[1].minimizer[j]) for j in 1:3] 
append!(MLEGGGH, OPTGGGH[1].minimizer[4:end]) 

MLEGGH = [exp(OPTGGH[1].minimizer[j]) for j in 1:2] 
append!(MLEGGH, OPTGGH[1].minimizer[3:end]) 

MLELNGH = [OPTLNGH[1].minimizer[1], exp(OPTLNGH[1].minimizer[2]), OPTLNGH[1].minimizer[3:end]...]

MLELLGH = [OPTLLGH[1].minimizer[1], exp(OPTLLGH[1].minimizer[2]), OPTLLGH[1].minimizer[3:end]...]

MLES = hcat(MLEPGWGH, MLEEWGH, MLEGGGH, [MLEGGH[1], MLEGGH[2], nothing, MLEGGH[3:end]...],
            [MLELNGH[1], MLELNGH[2], nothing, MLELNGH[3:end]...],[MLELLGH[1], MLELLGH[2], nothing, MLELLGH[3:end]...])

MLES = DataFrame(MLES, :auto)
 
rename!( MLES, ["PGWGH", "EWGH", "GGGH", "GGH", "LNGH", "LLGH"] )

# MLEs for GH models
println(MLES)

```

## Model Comparison

```@example 1
# AIC for models with PGW baseline hazard
AICPGWGH = 2*OPTPGWGH[1].minimum + 2*length(OPTPGWGH[1].minimizer)
AICPGWAFT = 2*OPTPGWAFT[1].minimum + 2*length(OPTPGWAFT[1].minimizer)
AICPGWPH = 2*OPTPGWPH[1].minimum + 2*length(OPTPGWPH[1].minimizer)
AICPGWAH = 2*OPTPGWAH[1].minimum + 2*length(OPTPGWAH[1].minimizer)

# AICs for models with LL baseline hazard
AICLLGH = 2*OPTLLGH[1].minimum + 2*length(OPTLLGH[1].minimizer)
AICLLAFT = 2*OPTLLAFT[1].minimum + 2*length(OPTLLAFT[1].minimizer)
AICLLPH = 2*OPTLLPH[1].minimum + 2*length(OPTLLPH[1].minimizer)
AICLLAH = 2*OPTLLAH[1].minimum + 2*length(OPTLLAH[1].minimizer)

# AICs for GH models with GG, EW, LN, and G hazards
AICGGGH = 2*OPTGGGH[1].minimum + 2*length(OPTGGGH[1].minimizer)
AICEWGH = 2*OPTEWGH[1].minimum + 2*length(OPTEWGH[1].minimizer)
AICLNGH = 2*OPTLNGH[1].minimum + 2*length(OPTLNGH[1].minimizer)
AICGGH = 2*OPTGGH[1].minimum + 2*length(OPTGGH[1].minimizer)



# All AICs
AICs = [AICPGWGH, AICPGWAFT, AICPGWPH, AICPGWAH,
          AICLLGH, AICLLAFT, AICLLPH, AICLLAH,
          AICGGGH, AICEWGH, AICLNGH, AICGGH]

println(AICs)

# Best model: LLGH
argmin(AICs)
```

## Baseline hazards for GH models

```@example 1
# Fitted baseline hazard functions for GH models
function PGWGHhaz(t::Float64) 
    out = hPGW(t, MLEPGWGH[1], MLEPGWGH[2], MLEPGWGH[3])  
    return out
end

function EWGHhaz(t::Float64) 
    out = hEW(t, MLEEWGH[1], MLEEWGH[2], MLEEWGH[3])  
    return out
end

function GGGHhaz(t::Float64) 
    out = hGenGamma(t, MLEGGGH[1], MLEGGGH[2], MLEGGGH[3])  
    return out
end

function GGHhaz(t::Float64) 
    out = hGamma(t, MLEGGH[1], MLEGGH[2])  
    return out
end

function LNGHhaz(t::Float64) 
    out = hLogNormal(t, MLELNGH[1], MLELNGH[2])  
    return out
end

function LLGHhaz(t::Float64) 
    out = hLogNormal(t, MLELLGH[1], MLELLGH[2])  
    return out
end

# Note that the baseline hazards associated to the top models look similar
plot(t -> PGWGHhaz(t),
      xlabel = "Time (years)", ylabel = "Baseline Hazard", title = "",
    xlims = (0.0001,maximum(times)),   xticks = 0:1:maximum(times), label = "", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue", ylims = (0,4))

plot!(t -> EWGHhaz(t),
   xlabel = "Time (years)", ylabel = "Baseline Hazard", title = "",
  xlims = (0.0001,maximum(times)),   xticks = 0:1:maximum(times), label = "", 
  xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
  xguidefontsize=18, yguidefontsize=18, linewidth=3,
  linecolor = "blue", ylims = (0,4), linestyle=:dash)

plot!(t -> GGGHhaz(t),
   xlabel = "Time (years)", ylabel = "Baseline Hazard", title = "",
  xlims = (0.0001,maximum(times)),   xticks = 0:1:maximum(times), label = "", 
  xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
  xguidefontsize=18, yguidefontsize=18, linewidth=3,
  linecolor = "blue", ylims = (0,4), linestyle=:dot)

plot!(t -> GGHhaz(t),
   xlabel = "Time (years)", ylabel = "Baseline Hazard", title = "",
  xlims = (0.0001,maximum(times)),   xticks = 0:1:maximum(times), label = "", 
  xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
  xguidefontsize=18, yguidefontsize=18, linewidth=3,
  linecolor = "blue", ylims = (0,4), linestyle=:dashdot)

plot!(t -> LNGHhaz(t),
      xlabel = "Time (years)", ylabel = "Baseline Hazard", title = "",
    xlims = (0.0001,maximum(times)),   xticks = 0:1:maximum(times), label = "", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue", ylims = (0,4), linestyle=:dashdotdot)

plot!(t -> LLGHhaz(t),
    xlabel = "Time (years)", ylabel = "Baseline Hazard", title = "",
  xlims = (0.0001,maximum(times)),   xticks = 0:1:maximum(times), label = "", 
  xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
  xguidefontsize=18, yguidefontsize=18, linewidth=3,
  linecolor = "red", ylims = (0,4), linestyle=:solid)
```

## Best-model summaries

```@example 1
# MLE in the original parameterisation
MLE = MLELLGH

println(MLE)

# 95% Confidence intervals under the reparameterisation
CI = ConfInt(FUN = OPTLLGH[2], MLE = OPTLLGH[1].minimizer, level = 0.95)

CI = DataFrame(CI, :auto)
 
rename!( CI, ["Lower", "Upper"] )

println(CI)
```

```@example 1
# Fitted baseline hazard function
plot(t -> LLGHhaz(t),
    xlabel = "Time (years)", ylabel = "Baseline Hazard", title = "Best Model",
  xlims = (0.0001,maximum(times)),   xticks = 0:1:maximum(times), label = "", 
  xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
  xguidefontsize=18, yguidefontsize=18, linewidth=3,
  linecolor = "red", ylims = (0,4), linestyle=:solid)

# Average population survival function and KM estimator

function pop_surv(t::Float64)
  p0 = size(des_t)[2]
  p1 = size(des)[2]
  theta1 = MLE[1]
  theta2 = MLE[2]
  alpha = MLE[3:(2+p0)]
  beta = MLE[(3+p0):(2+p0+p1)]
  x_alpha = des_t * alpha
  x_dif = des * beta - x_alpha
  out = mean( exp.( - chLogLogistic(t*exp.(x_alpha), theta1, theta2).* exp.(x_dif)  )  )
  return out
end


# Kaplan-Meier estimator 
km_fit = fit(KaplanMeier, times, df.cens)

ktimes = sort(unique(times))
ksurvival_probs = km_fit.survival

# Comparison
plot(ktimes, ksurvival_probs,
    xlabel = "Time (years)", ylabel = "Population Survival", title = "Best Model",
  xlims = (0.0001,maximum(times)),   xticks = 0:1:maximum(times), label = "", 
  xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
  xguidefontsize=18, yguidefontsize=18, linewidth=3,
  linecolor = "gray", ylims = (0,1), linestyle=:solid)

plot!(t -> pop_surv(t),
  xlabel = "Time (years)", ylabel = "Population Survival", title = "Best Model",
xlims = (0.0001,maximum(times)),   xticks = 0:1:maximum(times), label = "", 
xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
xguidefontsize=18, yguidefontsize=18, linewidth=3,
linecolor = "black", ylims = (0,1), linestyle=:solid)
```

```@example 1


# Confidence intervals for the survival function based on a normal approximation
# at specific time points t0

# Hessian and asymptotic covariance matrix
HESS = ForwardDiff.hessian(OPTLLGH[2], OPTLLGH[1].minimizer);
Sigma = inv(HESS);

#= 
A "hackish" workaround to a bug in ForwardDiff, which may
produce non-symmetric hessian matrices.
Here, I am replacing the lower diagonal of Sigma by its 
upper diagonal
=#
ps = size(Sigma, 1)
for i in 1:ps-1
    for j in i+1:ps
        Sigma[i, j] = Sigma[j, i]
    end
end


# Reparameterised MLE 
r_MLE = OPTLLGH[1].minimizer;


# The function to obtain approximate CIs based on Monte Carlo simulations 
# from the asymptotic normal distribution of the MLEs
# t0 : time where the confidence interval will be calculated
# level : confidence level
# nmc : number of Monte Carlo iterations

function ConfIntSurv(t0::Float64, level::Float64, nmc::Int64)
    p0 = size(des_t)[2]
    p1 = size(des)[2]
    mc = fill(0.0, nmc)
    function S_par(par::Vector{Float64})
        
        outs = mean( exp.( - chLogLogistic(t0.*exp.(des_t * par[3:(2+p0)]), par[1], par[2]) .*
                                        exp.( des * par[(3+p0):(2+p0+p1)].- des_t * par[3:(2+p0)])  )  )
                                        return outs
    end
    
    for i in 1:nmc
      mv_normal = MvNormal(r_MLE, Sigma) 
      val = rand(mv_normal, 1)
      val1 = [val[1],exp(val[2]),val[3:end]...]
      mc[i] = S_par(val1)
    end
    
    L = quantile(mc,(1-level)*0.5)
    U = quantile(mc,(1+level)*0.5)
    
    M = S_par(MLE)
    
    return [L,M,U]
end
  


# times for CIs calculations
timesCI = [1.0,2.5,5,7.5,10,12.5];

CIS = zeros(length(timesCI), 4);

for k in 1:length(timesCI)
CIS[k,:] = vcat(timesCI[k], ConfIntSurv(timesCI[k],0.95,10000))
end

CIS = DataFrame(CIS, :auto);
 
rename!( CIS, ["year","lower","population survival","upper"] );

println(CIS)

```

```@bibliography
Pages = ["HazRegJulia.md"]
Canonical = false
```