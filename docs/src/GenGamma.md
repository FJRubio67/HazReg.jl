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
using HazReg
```

# Examples

## Random number generation

```@example 1
#= Fix the seed =#
Random.seed!(123)
#= True values of the parameters =#
D = GeneralizedGamma(1,3,2) # sigma, nu, gamma
#= Simulation =#
sim = rand(D,1000);
```

## Some plots

```@example 1
#= Histogram and probability density function =#
histogram(sim, normalize=:pdf, color=:gray, 
          bins = range(0, 4, length=30), label = "")
plot!(t -> pdf(D,t),
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

plot(x -> ecdfsim(x), 0, 5, label = "ECDF", linecolor = "gray", linewidth=3)
plot!(t -> cdf(D,t),
      xlabel = "x", ylabel = "CDF vs. ECDF", title = "GenGamma distribution",
      xlims = (0,5),   xticks = 0:1:5, label = "CDF", 
      xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
      xguidefontsize=18, yguidefontsize=18, linewidth=3,
      linecolor = "blue")
```

```@example 1
#= Hazard function =#
plot(t -> haz(GeneralizedGamma(0.5, 1.5, 0.75),t),
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