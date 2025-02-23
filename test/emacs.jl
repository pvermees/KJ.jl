if !(@isdefined rerun)
    using Revise, Pkg
    Pkg.activate("/home/pvermees/git/KJ.jl")
    Pkg.instantiate()
    Pkg.precompile()
    cd("/home/pvermees/git/KJ.jl/test")
end

rerun = true

if false
    using KJ, Test, CSV, Infiltrator, DataFrames, Statistics
    import Plots

    snum = 23
    
    myrun = load("/home/pvermees/Documents/Plasmatrace/NHM/240708_PV_Zircon_Maps.csv",
                 "/home/pvermees/Documents/Plasmatrace/NHM/240708_PV_Zircon.Iolite.csv";
                 instrument="Agilent")
    deleteat!(myrun, 31)
    method = "U-Pb";
    channels = Dict("d"=>"Pb207",
                    "D"=>"Pb206",
                    "P"=>"U238");
    standards = Dict("91500_zr" => "91500")
    glass = Dict("NIST612" => "NIST612")
    blk, fit = process!(myrun,method,channels,standards,glass,
                        nblank=2,ndrift=2,ndown=0)
    export2IsoplotR(myrun,method,channels,blk,fit;
                    fname="/home/pvermees/temp/NHM.json")
    p = plot(myrun[snum],method,channels,blk,fit,standards,glass;
             transformation="log",den=nothing)
else
    include("runtests.jl")
end
