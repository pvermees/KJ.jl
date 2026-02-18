using Documenter, KJ

makedocs(
    ;
    modules = [KJ],
    sitename="KJ.jl",
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    ;
    repo = "github.com/pvermees/KJ.jl",
    devbranch = "beta",
    # Only push previews if all the relevant environment variables are non-empty. This is an
    # attempt to work around https://github.com/JuliaDocs/Documenter.jl/issues/2048.
    push_preview = all(!isempty, (get(ENV, "GITHUB_TOKEN", ""), get(ENV, "DOCUMENTER_KEY", ""))),
)
