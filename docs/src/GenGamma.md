# Generalised Gamma Distribution

# The Generalised Gamma Distribution

The [Generalised Gamma](https://en.wikipedia.org/wiki/Generalized_gamma_distribution) (GG) distribution [stacy:1962](@cite) is a three-parameter distribution with support on ${\mathbb R}_+$. The corresponding hazard function can accommodate bathtub, unimodal and monotone (increasing and decreasing) hazard shapes. The GG distribution has become popular in survival analysis due to its flexibility. Other flexible distributions that can account for these hazard shapes are discussed in @rubio:2021 and @jones:2015. 

## Probability Density Function

The pdf of the GG distribution is

$$f(t;\theta,\kappa,\delta) = \dfrac{\delta}{\Gamma\left(\frac{\kappa}{\delta}\right)\theta^\kappa}  t^{{\kappa-1}}e^{{-\left(\frac{t}{\theta}\right)^{\delta}}},$$
where $\theta>0$ is a scale parameter, and $\kappa,\delta >0$ are shape parameters.

## Cumulative Distribution Function

The CDF of the GG distribution is

$$F(t;\theta,\kappa,\delta) = {\frac  {\gamma \left( \frac{\kappa}{\delta},\left(\frac{t}{\theta}\right)^{\delta}\right)}{\Gamma\left(\frac{\kappa}{\delta}\right)}},$$

where where $\gamma (\cdot )$ denotes the lower incomplete gamma function. The survival function can be obtained using the relationship $S(t;\theta,\kappa,\delta)=1-F(t;\theta,\kappa,\delta)$. An interesting relationship between the [Gamma CDF](https://en.wikipedia.org/wiki/Gamma_distribution) ($G(t;\theta,\kappa)$, scale $\theta$ and shape $\kappa$) and the GG CDF is

$$F(t;\theta,\kappa,\delta) = G\left(t^\delta; \theta^\delta, \frac{\kappa}{\delta}\right).$$
This allows the implementation of the GG CDF using the Julia command `Gamma`.

###  Hazard Function

The hazard function of the GG distribution is

$$h(t;\theta,\kappa,\delta) = \dfrac{f(t;\theta,\kappa,\delta)}{1-F(t;\theta,\kappa,\delta)}.$$

The survival function can be obtained as $S(t;\theta,\kappa,\delta)=1-F(t;\theta,\kappa,\delta)$, and the cumulative hazard function as $H(t;\theta,\kappa,\delta) = -\log S(t;\theta,\kappa,\delta)$, as usual. The connection of the GG CDF with the Gamma distribution allows for writing these functions in terms of the Julia command `Gamma` as shown in the following code.

The following Julia code shows the implementation of the pdf, survival function, hazard function, cumulative hazard function, quantile function, and random number generation associated to the Generalised Gamma distribution using the Julia package `HazReg.jl`. Some illustrative examples are also presented.


See also: 

- [HazReg.jl Julia Package](https://github.com/FJRubio67/HazReg.jl) 

- [The Generalised Gamma Distribution](https://rpubs.com/FJRubio/GG)

# Required packages
```@example 1
using Distributions
using Random
using Plots
using StatsBase
using SpecialFunctions
```

# Functions

```@example 1
#=
Generalised Gamma (GG) Distribution: Hazard, cumulative hazard,
probability density function, random number generation, and survival function
=#


#= Generalised Gamma (GG) probability density function.
t: positive argument

sigma: scale parameter

nu: shape parameter

gamma: shape parameter

logd: log scale (true or false)
=#

function pdfGenGamma(t, sigma, nu, gamma, logd::Bool = false)
    val = log(gamma) .- nu * log(sigma) .- loggamma(nu/gamma) .+ 
    (nu - 1) * log.(t) .- (t/sigma).^gamma
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

logp: log scale (true or false)
=#

function ccdfGenGamma(t, sigma, nu, gamma, logp::Bool = false)
    tp = t.^gamma
    val = logccdf.(Gamma(nu/gamma, sigma^gamma), tp)
    if logp
        return val
    else
        return exp.(val)
    end
end


#= Generalised Gamma (GG) hazard function.
t: positive argument

sigma: scale parameter

nu: shape parameter

gamma: shape parameter

logh: log scale (true or false)
=#

function hGenGamma(t, sigma, nu, gamma, logh::Bool = false)
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
    out = quantile.(Gamma(nu/gamma, sigma^gamma), p).^(1/gamma)
    return out
end
```

# Examples

## Random number generation

```@example 1
#= Fix the seed =#
Random.seed!(123)
#= True values of the parameters =#
sigma0 = 1
nu0 = 3
gamma0 = 2
#= Simulation =#
sim = randGenGamma(1000, sigma0, nu0, gamma0);
```

## Some plots

```@example 1
#= Histogram and probability density function =#
histogram(sim, normalize=:pdf, color=:gray, 
          bins = range(0, 4, length=30), label = "")
plot!(t -> pdfGenGamma(t, sigma0, nu0, gamma0),
      xlabel = "x", ylabel = "Density", title = "GenGamma distribution",
    xlims = (0,4),   xticks = 0:1:4, label = "", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue")
```


```@example 1
#= Empirical CDF and CDF =#

#= Empirical CDF=#
ecdfsim = ecdf(sim)

#= ad hoc CDF =#
function cdfGenGamma(t, sigma, nu, gamma) 
val = 1 .- ccdfGenGamma.(t, sigma, nu, gamma)
return val
end

plot(x -> ecdfsim(x), 0, 5, label = "ECDF", linecolor = "gray", linewidth=3)
plot!(t -> cdfGenGamma(t, sigma0, nu0, gamma0),
      xlabel = "x", ylabel = "CDF vs. ECDF", title = "GenGamma distribution",
      xlims = (0,5),   xticks = 0:1:5, label = "CDF", 
      xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
      xguidefontsize=18, yguidefontsize=18, linewidth=3,
      linecolor = "blue")
```

```@example 1
#= Hazard function =#
plot(t -> hGenGamma(t, 0.5, 1.5, 0.75),
     xlabel = "x", ylabel = "Hazard", title = "GenGamma distribution",
     xlims = (0,15),   xticks = 0:1:15, label = "", 
     xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
     xguidefontsize=18, yguidefontsize=18, linewidth=3,
     linecolor = "blue")
```


```@bibliography
Pages = ["GenGamma.md"]
Canonical = false
```