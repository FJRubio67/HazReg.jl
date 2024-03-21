module HazReg

#= Imports =#

using Distributions
using Random
using StatsBase
using Optim
using LinearAlgebra
using SpecialFunctions
using ForwardDiff
import Base.rand

#= Includes =#
include("utils.jl")
include("univ_dists/ExponentiatedWeibull.jl")
include("univ_dists/GeneralizedGamma.jl")
include("univ_dists/PowerGeneralizedWeibull.jl")
include("univ_dists/LogLogistic.jl")
include("GHMLE.jl")
include("simGH.jl")
include("ConfInt.jl")

#= Exports =#
export  GHMLE
export  simGH
export  ConfInt
export  standardise

export ExponentiatedWeibull
export GeneralizedGamma
export PowerGeneralizedWeibull
export LogLogistic
export haz, loghaz, cumhaz

end
