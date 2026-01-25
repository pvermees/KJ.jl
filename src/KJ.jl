"""
KJ

A Julia package to process LA-ICP-MS data
"""
module KJ

using Infiltrator, DataFrames, Dates, LinearAlgebra, Printf, Plots

import Statistics, Distributions, Optim, CSV, ForwardDiff, FiniteDiff

include("init.jl")
include("errors.jl")
include("json.jl")
include("types.jl")
include("OrderedDict.jl")
include("Pairing.jl")
include("Calibration.jl")
include("KJmethod.jl")
include("KJfit.jl")
include("io.jl")
include("accessors.jl")
include("toolbox.jl")
include("windows.jl")
include("plots.jl")
include("blank.jl")
include("outliers.jl")
include("MCD.jl")
include("interference.jl")
# include("blocks.jl")
# include("anchors.jl")
# include("parser.jl")
# include("crunch.jl")
# include("bias.jl")
# include("fractionation.jl")
# include("process.jl")
# include("atomic.jl")
# include("concentrations.jl")
# include("averat.jl")
# include("internochron.jl")
# include("internoplot.jl")
include("TUImessages.jl")
include("TUIactions.jl")
include("TUI.jl")

init_KJ!()

end
