# Some hazard functions in Julia

# Hazard and cumulative hazard functions

The hazard and the cumulative hazard functions play a crucial role in survival analysis. These functions define the likelihood function in the presence of censored observations. Thus, they are important in many context. For more information about these functions, see [Short course on Parametric Survival Analysis
](https://github.com/FJRubio67/ShortCourseParamSurvival).

This tutorial contains the implementation of the hazard and cumulative hazard functions for some popular distributions using the Julia package `HazReg.jl`, such as the LogNormal, LogLogistic, Weibull, and Gamma distributions. 

The Julia implementation of other, less common, distributions that are of interest in survival analysis can be found in the Julia package `HazReg.jl`: [`HazReg.jl`](https://github.com/FJRubio67/HazReg.jl).

See also: 

- [HazReg.jl Julia Package](https://github.com/FJRubio67/HazReg.jl) 

# Required packages
```@example 1
using Distributions
using Plots
using StatsBase
using HazReg
```

# Examples

## LogNormal

```@example 1
#= Hazard function =#
plot(t -> hLogNormal(t, 0.5, 1),
      xlabel = "x", ylabel = "Hazard", title = "LogNormal distribution",
    xlims = (0,10),   xticks = 0:1:10, label = "", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue")
```

```@example 1
#= Cumulative Hazard function =#
plot(t -> chLogNormal(t, 0.5, 1),
      xlabel = "x", ylabel = "Cumulative Hazard", title = "LogNormal distribution",
    xlims = (0,10),   xticks = 0:1:10, label = "", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue")
```


## LogLogistic

```@example 1
#= Hazard function =#
plot(t -> hLogLogistic(t, 1, 0.5),
      xlabel = "x", ylabel = "Hazard", title = "LogLogistic distribution",
    xlims = (0,10),   xticks = 0:1:10, label = "", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue")
```

```@example 1
#= Cumulative Hazard function =#
plot(t -> chLogLogistic(t, 1, 0.5),
      xlabel = "x", ylabel = "Cumulative Hazard", title = "LogLogistic distribution",
    xlims = (0,10),   xticks = 0:1:10, label = "", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue")
```

## Weibull

```@example 1
#= Weibull function =#
plot(t -> hWeibull(t, 3, 0.5),
      xlabel = "x", ylabel = "Hazard", title = "Weibull distribution",
    xlims = (0,10),   xticks = 0:1:10, label = "", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue")
```

```@example 1
#= Cumulative Hazard function =#
plot(t -> chWeibull(t, 3, 0.5),
      xlabel = "x", ylabel = "Cumulative Hazard", title = "Weibull distribution",
    xlims = (0,10),   xticks = 0:1:10, label = "", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue")
```


## Gamma

```@example 1
#= Weibull function =#
plot(t -> hGamma(t, 3, 0.5),
      xlabel = "x", ylabel = "Hazard", title = "Gamma distribution",
    xlims = (0,10),   xticks = 0:1:10, label = "", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue")
```

```@example 1
#= Cumulative Hazard function =#
plot(t -> chGamma(t, 3, 0.5),
      xlabel = "x", ylabel = "Cumulative Hazard", title = "Gamma distribution",
    xlims = (0,10),   xticks = 0:1:10, label = "", 
    xtickfont = font(16, "Courier"),  ytickfont = font(16, "Courier"),
    xguidefontsize=18, yguidefontsize=18, linewidth=3,
    linecolor = "blue")
```