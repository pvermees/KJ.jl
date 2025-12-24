"""
KJ

A Julia package to process LA-ICP-MS data
"""
module KJ

using Infiltrator, DataFrames, Dates, LinearAlgebra, Printf

import Plots, Statistics, Distributions, Optim, CSV, ForwardDiff, FiniteDiff

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
include("interference.jl")
include("fractionation.jl")
include("process.jl")
include("atomic.jl")
include("concentrations.jl")
include("averat.jl")
include("plots.jl")
include("internochron.jl")
include("internoplot.jl")
include("TUImessages.jl")
include("TUIactions.jl")
include("TUI.jl")

init_KJ!()

end
