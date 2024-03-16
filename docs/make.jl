using Documenter
using HazReg

DocMeta.setdocmeta!(HazReg, :DocTestSetup, :(using HazReg); recursive=true)

makedocs(
    sitename = "HazReg.jl",
    format = Documenter.HTML(),
    modules = [HazReg]
)

deploydocs(;
    repo="github.com/FJRubio67/HazReg.jl",
    devbranch="main",
)
