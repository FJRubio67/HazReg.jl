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
include("univ_dists/ExponentiatedWeibull.jl")
include("univ_dists/GeneralizedGamma.jl")
include("univ_dists/PowerGeneralizedWeibull.jl")
include("univ_dists/other_hazards.jl")
include("GHMLE.jl")
include("simGH.jl")
include("ConfInt.jl")
include("standardise.jl")

#= Exports =#
export  GHMLE
export  simGH
export  ConfInt
export  standardise

export ExponentiatedWeibull
export GeneralizedGamma
export PowerGeneralizedWeibull
export haz, loghaz, cumhaz

end
