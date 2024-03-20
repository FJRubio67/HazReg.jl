using HazReg
using Test

@testset "HazReg.jl" begin
   haz(LogNormal(0.0, 1.0), 1.0)
end
