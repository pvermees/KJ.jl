module KJ

using Dates, DataFrames, Printf, Infiltrator, LinearAlgebra
import Plots, Statistics, Optim, CSV

include("errors.jl")
include("json.jl")
include("types.jl")
include("accessors.jl")
include("toolbox.jl")
include("io.jl")
include("plots.jl")
include("crunch.jl")
include("fractionation.jl")
include("process.jl")
include("averat.jl")
include("internochron.jl")
include("TUImessages.jl")
include("TUIactions.jl")
include("TUI.jl")

init_KJ!()

end
