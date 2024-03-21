using HazReg
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(HazReg, :DocTestSetup, :(using HazReg); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__,"src","assets","references.bib"),
    style=:numeric
)

makedocs(;
    plugins=[bib],
    modules = [HazReg],
    authors="F.J.Rubio",
    repo="https://github.com/FJRubio67/HazReg.jl/blob/{commit}{path}#{line}",
    sitename = "HazReg.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://FJRubio67.github.io/HazReg.jl",
        assets=String["assets/citations.css"],
        collapselevel=3,
    ),
)

deploydocs(;
    repo="github.com/FJRubio67/HazReg.jl",
    devbranch="main",
)
