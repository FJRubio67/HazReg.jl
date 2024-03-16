# Power Generalised Weibull Distribution

# The Power Generalised Weibull Distribution

The Power Generalised Weibull (PGW) distribution [nikulin:2009](@cite) is a three-parameter distribution with support on ${\mathbb R}_+$. The corresponding hazard function can accommodate bathtub, unimodal and monotone (increasing and decreasing) hazard shapes. The PGW distribution has become popular in survival analysis given the tractability of its hazard and survival functions. Other flexible distributions that can account for these hazard shapes are discussed in @rubio:2021 and @jones:2015. 

## Probability Density Function

The pdf of the PGW distribution is

$$f(t;\sigma,\nu,\gamma) = \dfrac{\nu}{\gamma \sigma^\nu}t^{\nu-1} \left[ 1 + \left(\dfrac{t}{\sigma}\right)^\nu\right]^{\left(\frac{1}{\gamma}-1\right)} \exp\left\{ 1- \left[ 1 + \left(\dfrac{t}{\sigma}\right)^\nu\right]^{\frac{1}{\gamma}}
\right\},$$

where $\sigma>0$ is a scale parameter, and $\nu,\gamma >0$ are shape parameters.

## Survival Function

The survival function of the PGW distribution is

$$S(t;\sigma,\nu,\gamma) = \exp\left\{ 1- \left[ 1 + \left(\dfrac{t}{\sigma}\right)^\nu\right]^{\frac{1}{\gamma}}
\right\}.$$

##  Hazard Function

The hazard function of the PGW distribution is

$$h(t;\sigma,\nu,\gamma) = \dfrac{\nu}{\gamma \sigma^\nu}t^{\nu-1} \left[ 1 + \left(\dfrac{t}{\sigma}\right)^\nu\right]^{\left(\frac{1}{\gamma}-1\right)}.$$
The cdf can be obtained as $F(t;\sigma,\nu,\gamma)=1-S(t;\sigma,\nu,\gamma)$, and the cumulative hazard function as $H(t;\sigma,\nu,\gamma) = -\log S(t;\sigma,\nu,\gamma)$, as usual.

## Quantile Function

The quantile function of the PGW distribution is

$$Q(p;\sigma,\nu,\gamma) = \sigma \left[ \left( 1 - \log(1-p) \right)^{\gamma} - 1 \right]^{\frac{1}{\nu}},$$

where $p\in(0,1)$.

The following Julia code shows the implementation of the pdf, survival function, hazard function, cumulative hazard function, quantile function, and random number generation associated to the PGW distribution using the Julia package `HazReg.jl`. Some illustrative examples are also presented.


See also: 

- [HazReg.jl Julia Package](https://github.com/FJRubio67/HazReg.jl) 

- [The Power Generalised Weibull Distribution](https://rpubs.com/FJRubio/PGW)

# Required packages
```@example 1
using Distributions
using Random
using Plots
using StatsBase
using HazReg
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
sim = randPGW(1000, sigma0, nu0, gamma0);
```

## Some plots

```@example 1
#= Histogram and probability density function =#
histogram(sim, normalize=:pdf, color=:gray, 
          bins = range(0, 5, length=30), label = "")
plot!(t -> pdfPGW(t, sigma0, nu0, gamma0),
      xlabel = "x", ylabel = "Density", title = "PGW distribution",
    xlims = (0,5),   xticks = 0:1:5, label = "", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue")
```

```@example 1
#= Empirical CDF and CDF =#

#= Empirical CDF=#
ecdfsim = ecdf(sim)

#= ad hoc CDF =#
function cdfPGW(t, sigma, nu, gamma) 
        val = 1 .- ccdfPGW.(t, sigma, nu, gamma)
        return val
end

plot(x -> ecdfsim(x), 0, 5, label = "ECDF", linecolor = "gray", linewidth=3)
plot!(t -> cdfPGW(t, sigma0, nu0, gamma0),
      xlabel = "x", ylabel = "CDF vs. ECDF", title = "PGW distribution",
    xlims = (0,5),   xticks = 0:1:5, label = "CDF", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue")
```


```@example 1
#= Hazard function =#
plot(t -> hPGW(t, 0.5, 2, 5),
      xlabel = "x", ylabel = "Hazard", title = "PGW distribution",
    xlims = (0,10),   xticks = 0:1:10, label = "", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue")
```



```@bibliography
Pages = ["PGW.md"]
Canonical = false
```