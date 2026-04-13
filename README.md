# HazReg.jl

## Parametric Hazard-Based Regression Models for Survival Data

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://FJRubio67.github.io/HazReg.jl/dev)
[![Julia package](https://img.shields.io/badge/language-Julia-9558B2.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-v0.1.0-green.svg)](https://github.com/FJRubio67/HazReg.jl/releases/tag/v0.1.0%2Bdoc1)

## Overview

`HazReg.jl` is a Julia package for fitting parametric hazard-based regression
models for overall survival data. It is the Julia counterpart of the
[HazReg R package](https://github.com/FJRubio67/HazReg), and is built around
the **General Hazard (GH)** structure, which nests the most widely used hazard
regression models as special cases:

| Model | Abbreviation | Special case of GH |
|---|---|---|
| General Hazard | GH | — |
| Proportional Hazards | PH | ✓ |
| Accelerated Failure Time | AFT | ✓ |
| Accelerated Hazards | AH | ✓ |

Models are fitted by maximum likelihood using the
[`Optim.jl`](https://github.com/JuliaNLSolvers/Optim.jl) package. Users should
specify initial values and verify convergence, as is standard practice. The
following optimisation methods are available:

| Method key | Algorithm |
|---|---|
| `"NM"` | Nelder-Mead |
| `"N"` | Newton |
| `"LBFGS"` | L-BFGS |
| `"CG"` | Conjugate Gradient |
| `"GD"` | Gradient Descent |

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/FJRubio67/HazReg.jl")
using HazReg
```

## Main function

| Function | Description |
|---|---|
| `GHMLE` | Maximum likelihood estimation for GH, PH, AFT, and AH models |
| `simGH` | Simulation of survival times from a GH structure |

Full documentation is available at the
[dev docs](https://FJRubio67.github.io/HazReg.jl/dev).

## Supported baseline hazard distributions

All positive parameters are log-transformed for unconstrained optimisation.

| Distribution | Key | GH | PH | AFT | AH |
|---|---|:---:|:---:|:---:|:---:|
| [Power Generalised Weibull](http://rpubs.com/FJRubio/PGW) | PGW | ✓ | ✓ | ✓ | ✓ |
| [Exponentiated Weibull](http://rpubs.com/FJRubio/EWD) | EW | ✓ | ✓ | ✓ | ✓ |
| [Generalised Gamma](http://rpubs.com/FJRubio/GG) | GenGamma | ✓ | ✓ | ✓ | ✓ |
| [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution) | Gamma | ✓ | ✓ | ✓ | ✓ |
| [Log-normal](https://en.wikipedia.org/wiki/Log-normal_distribution) | LogNormal | ✓ | ✓ | ✓ | ✓ |
| [Log-logistic](https://en.wikipedia.org/wiki/Log-logistic_distribution) | LogLogistic | ✓ | ✓ | ✓ | ✓ |
| [Weibull](https://en.wikipedia.org/wiki/Weibull_distribution) | Weibull | — | ✓ | ✓ | ✓ |

## Tutorials and examples

- [HazReg.jl: Parametric hazard-based regression models](https://fjrubio.quarto.pub/hazregjulia/) — main illustrative tutorial
- [Simulating survival times from a GH structure](https://fjrubio.quarto.pub/simghjulia/) — simulation guide
- [Power Generalised Weibull distribution](https://fjrubio.quarto.pub/power-generalised-weibull-distribution/) — distribution reference
- [Generalised Gamma distribution](https://fjrubio.quarto.pub/generalised-gamma-distribution/) — distribution reference
- [Exponentiated Weibull distribution](https://fjrubio.quarto.pub/exponentiated-weibull-distribution/) — distribution reference
- [Some hazard functions in Julia](https://fjrubio.quarto.pub/some-hazard-functions-in-julia/) — companion notes

## Related resources

- [HazReg](https://github.com/FJRubio67/HazReg) — R implementation of the same models
- [PTCMGH](https://github.com/FJRubio67/PTCMGH) — Promotion Time Cure Models with a GH structure (R)
- [SurvMODE](https://github.com/FJRubio67/SurvMODE) — survival modelling via ODEs (R and Julia)
- [Short course on Parametric Survival Analysis](https://github.com/FJRubio67/ShortCourseParamSurvival)

## License

This package is licensed under the [MIT License](LICENSE).
