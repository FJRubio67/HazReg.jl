# Exponentiated Weibull Distribution

## The Exponentiated Weibull distribution

The [Exponentiated Weibull distribution](https://en.wikipedia.org/wiki/Exponentiated_Weibull_distribution) is a generalisation of the [Weibull distribution](https://en.wikipedia.org/wiki/Weibull_distribution) which is obtained by exponentiating the Weibull cumulative distribution function. This simple transformation adds a second shape parameter that, interestingly, induces a lot of flexibility on the hazard function. The hazard function of the Exponentiated Weibull distribution can capture the basic shapes: constant, increasing, decreasing, bathtub, and unimodal, making it appealing for survival models.

The probability density function and cumulative distribution function of the Exponentiated Weibull distribution are respectively given by:

$$\begin{split}
f_{EW}(t) &=  \alpha \dfrac{\kappa}{\lambda} \left(\dfrac{t}{\lambda}\right)^{\kappa-1} \left[1-\exp\left\{-\left(\dfrac{t}{\lambda}\right)^{\kappa}\right\}\right]^{\alpha-1} \exp\left\{-\left(\dfrac{t}{\lambda}\right)^{\kappa}\right\}, \\
F_{EW}(t) &= \left[1-\exp\left\{-\left(\dfrac{t}{\lambda}\right)^{\kappa}\right\}\right]^{\alpha},
\end{split}$$

where $t>0$, $\alpha>0$, $\lambda>0$, and $\kappa>0$. The parameter $\lambda$ is a scale parameter, $\kappa$ is a shape parameter, and $\alpha$ is the power (shape) parameter. 

The following Julia code contains an implementation of the probability density function, cumulative distribution function, quantile function, random number generation, hazard function, and cumulative hazard function of the Exponentiated Weibull distribution using the corresponding Julia package `HazReg.jl`.
Some illustrative examples are also presented.

See also: 

- [HazReg.jl Julia Package](https://github.com/FJRubio67/HazReg.jl) 

- [The Exponentiated Weibull distribution](https://rpubs.com/FJRubio/EWD)

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
sim = randEW(1000, sigma0, nu0, gamma0);
```

## Some plots

```@example 1
#= Histogram and probability density function =#
histogram(sim, normalize=:pdf, color=:gray, 
          bins = range(0, 3, length=30), label = "")
plot!(t -> pdfEW(t, sigma0, nu0, gamma0),
      xlabel = "x", ylabel = "Density", title = "EW distribution",
    xlims = (0,3),   xticks = 0:1:3, label = "", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue")
```


```@example 1
#= Empirical CDF and CDF =#

#= Empirical CDF=#
ecdfsim = ecdf(sim)

plot(x -> ecdfsim(x), 0, 3, label = "ECDF", linecolor = "gray", linewidth=3)
plot!(t -> cdfEW(t, sigma0, nu0, gamma0),
      xlabel = "x", ylabel = "CDF vs. ECDF", title = "EW distribution",
      xlims = (0,3),   xticks = 0:1:3, label = "CDF", 
      xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
      xguidefontsize=18, yguidefontsize=18, linewidth=3,
      linecolor = "blue")
```

```@example 1
#= Hazard function =#
plot(t -> hEW(t, 0.25, 0.5, 5),
     xlabel = "x", ylabel = "Hazard", title = "EW distribution",
     xlims = (0,10),   xticks = 0:1:10, label = "", 
     xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
     xguidefontsize=18, yguidefontsize=18, linewidth=3,
     linecolor = "blue")
```