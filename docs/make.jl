using Documenter
using HazReg

makedocs(
    sitename = "HazReg",
    format = Documenter.HTML(),
    modules = [HazReg]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
