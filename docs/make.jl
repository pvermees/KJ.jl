using Documenter, KJ

makedocs(
    ;
    modules = [KJ],
    sitename="KJ.jl",
    pages = [
        "Home" => "index.md",
    ],
)
