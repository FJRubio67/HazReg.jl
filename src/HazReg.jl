module HazReg

#= Imports =#

using Distributions
using Random
using StatsBase
using Optim
using LinearAlgebra
using SpecialFunctions
using ForwardDiff

#= Includes =#
include("hazards.jl")
include("EW.jl")
include("GenGamma.jl")
include("PGW.jl")
include("GHMLE.jl")
include("simGH.jl")
include("ConfInt.jl")
include("standardise.jl")

#= Exports =#
export  GHMLE
export  simGH
export  ConfInt
export  standardise

export  pdfEW
export  cdfEW
export  qEW
export  hEW
export  chEW
export  randEW

export  pdfPGW
export  ccdfPGW
export  qPGW
export  hPGW
export  chPGW
export  randPGW

export  ccdfGenGamma
export  pdfGenGamma
export  qGenGamma
export  hGenGamma
export  chGenGamma
export  randGenGamma

export hLogNormal
export chLogNormal

export hLogLogistic
export chLogLogistic

export hWeibull
export chWeibull

export hGamma
export chGamma

end
