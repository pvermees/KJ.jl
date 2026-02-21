using Documenter, KJ

makedocs(
    ;
    modules = [KJ],
    sitename = "KJ.jl",
    pages = [
        "Home" => "index.md",
        "TUI" => "tui.md",
        "REPL" => "repl.md",
        "Developers" => "developers.md",
        "API" => "api.md",
    ],
    # Ensure source links resolve even if the module is not loaded from a git checkout.
    remotes = Dict(
        joinpath(@__DIR__, "..") => Documenter.Remotes.GitHub("pvermees", "KJ.jl"),
    ),
    format = Documenter.HTML(
        edit_link = "main",
        collapselevel = 3,
        ),
)

deploydocs(
    ;
    repo = "github.com/pvermees/KJ.jl.git",
    devbranch = "main",
    push_preview = false,
)
