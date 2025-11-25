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
include("KJmethod.jl")
include("KJfit.jl")
include("accessors.jl")
include("outliers.jl")
include("MCD.jl")
include("blocks.jl")
include("io.jl")
include("anchors.jl")
include("windows.jl")
include("toolbox.jl")
include("init.jl")
include("parser.jl")
include("crunch.jl")
include("blank.jl")
include("fractionation.jl")
include("process.jl")
include("averat.jl")
include("plots.jl")
include("internochron.jl")
include("internoplot.jl")
include("TUImessages.jl")
include("TUIactions.jl")
include("TUI.jl")

init_KJ!()

end
