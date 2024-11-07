if !(@isdefined rerun)
    using Revise, Pkg
    Pkg.activate("/home/pvermees/git/Plasmatrace.jl")
    Pkg.instantiate()
    Pkg.precompile()
    cd("/home/pvermees/git/Plasmatrace.jl/test")
end

rerun = true

include("runtests.jl")
