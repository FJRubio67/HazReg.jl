using HazReg
using Test
using Distributions

@testset "HazReg.jl" begin
   haz(LogNormal(0.0, 1.0), 1.0)
end
