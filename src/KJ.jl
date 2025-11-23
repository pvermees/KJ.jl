"""
KJ

A Julia package to process LA-ICP-MS data
"""
module KJ

using Dates, DataFrames, Printf, Infiltrator, LinearAlgebra, ForwardDiff
import Plots, Statistics, Distributions, Optim, CSV

include("errors.jl")
include("json.jl")
include("types.jl")
include("OrderedDict.jl")
include("accessors.jl")
include("windows.jl")
include("blocks.jl")
include("toolbox.jl")
include("init.jl")
include("io.jl")
include("parser.jl")
include("plots.jl")
include("crunch.jl")
include("outliers.jl")
include("MCD.jl")
include("blank.jl")
include("fractionation.jl")
include("process.jl")
include("averat.jl")
include("internochron.jl")
include("internoplot.jl")
include("TUImessages.jl")
include("TUIactions.jl")
include("TUI.jl")

init_KJ!()

end
