using Documenter
using KJ

makedocs(
    sitename = "KJ",
    format = Documenter.HTML(),
    modules = [KJ],
    pages    = [
        "Home" => "index.md",
        "API"  => "api.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
