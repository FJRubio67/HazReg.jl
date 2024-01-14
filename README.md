# HazReg.jl

The `HazReg` Julia package implements the following parametric hazard-based regression models for (overall) survival data.

- General Hazard (GH) model.

- Accelerated Failure Time (AFT) model.

- Proportional Hazards (PH) model.

- Accelerated Hazards (AH) model.


These models are fitted using the Julia package `Optim`, with methods "NM" (NelderMead), "N" (Newton), "LBFGS" (LBFGS), "CG" (ConjugateGradient), "GD" (GradientDescent). Thus, the user needs to specify the initial points and to check the convergence of the optimisation step, as usual.

To install the package, use:

```
Pkg.add(url="https://github.com/FJRubio67/HazReg.jl")

using HazReg
```


The current version of the `HazReg` Julia package implements the following parametric baseline hazards for the models discussed in the previous section, using the command `GHMLE`.

- [Power Generalised Weibull](http://rpubs.com/FJRubio/PGW) (PGW) distribution. 
 
- [Exponentiated Weibull](http://rpubs.com/FJRubio/EWD) (EW) distribution. 
 
- [Generalised Gamma](http://rpubs.com/FJRubio/GG) (GenGamma) distribuiton. 

- [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution) (Gamma) distribution. 

- [Lognormal](https://en.wikipedia.org/wiki/Log-normal_distribution) (LogNormal) distribution. 

- [Log-logistic](https://en.wikipedia.org/wiki/Log-logistic_distribution) (LogLogistic) distribution. 

- [Weibull](https://en.wikipedia.org/wiki/Weibull_distribution) (Weibull) distribution. (only for AFT, PH, and AH models) 


All positive parameters are transformed into the real line using a `log` link (reparameterisation).

Illustrative examples and a description of the available models can be found at:

1. [HazReg.jl: Parametric Hazard-based regression models for survival data](https://fjrubio.quarto.pub/hazregjulia/)
2. [simGHJulia](https://fjrubio.quarto.pub/simghjulia/)
3. [Power Generalised Weibull Distribution](https://fjrubio.quarto.pub/power-generalised-weibull-distribution/)
4. [Generalised Gamma Distribution](https://fjrubio.quarto.pub/generalised-gamma-distribution/)
5. [Exponentiated Weibull Distribution](https://fjrubio.quarto.pub/exponentiated-weibull-distribution/)
6. [Some hazard functions in Julia](https://fjrubio.quarto.pub/some-hazard-functions-in-julia/)


See also: 
- [HazReg R package](https://github.com/FJRubio67/HazReg)

