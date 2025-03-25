if !(@isdefined rerun)
    using Revise, Pkg
    Pkg.activate("/home/pvermees/git/KJ.jl")
    Pkg.instantiate()
    Pkg.precompile()
    cd("/home/pvermees/git/KJ.jl/test")
end

rerun = true

option = "NHM-carbonate" # "runtests" # "KJgui" # "Abdulkadir" # "NHM" #

if option != "runtests"
    using KJ, Test, CSV, Infiltrator, DataFrames, Statistics, Plots
end

if option == "Abdulkadir"

    myrun = load("/home/pvermees/Dropbox/Plasmatrace/Abdulkadir";
                 instrument="Agilent")

    snames = summarise(myrun).name
    snums = findall(x -> occursin("MAD", x), snames)

    method = "U-Pb"
    channels = Dict("d"=>"Pb207",
                    "D"=>"Pb206",
                    "P"=>"U238")
    standards = Dict("MAD_ap" => "MAD")
    glass = Dict("NIST612" => "GLASS")
    sett0!(myrun,6.5)
    setBwin!(myrun)
    setSwin!(myrun)
    blk, fit = process!(myrun,method,channels,standards,glass,
                        nblank=1,ndrift=1,ndown=1)
    p = KJ.plot(myrun[166],method,channels,blk,fit,standards,glass;
                transformation="log")
    display(p)
    if true
        export2IsoplotR(myrun,method,channels,blk,fit;
                        fname="/home/pvermees/temp/Abdulkadir.json")
        rm("/home/pvermees/temp/Abdulkadir.pdf")
        for snum in snums
            savefig(KJ.plot(myrun[snum],method,channels,blk,fit,standards,glass;
                            transformation="log"),#,den="Pb206"),
                    "/home/pvermees/temp/temp.pdf")
            append_pdf!("/home/pvermees/temp/Abdulkadir.pdf",
                        "/home/pvermees/temp/temp.pdf",
                        cleanup=true)
        end
    end
elseif option == "NHM"
    snum = 23
    myrun = load("/home/pvermees/Documents/Plasmatrace/NHM/240708_PV_Zircon_Maps.csv",
                 "/home/pvermees/Documents/Plasmatrace/NHM/240708_PV_Zircon.Iolite.csv";
                 instrument="Agilent")
    deleteat!(myrun, 31)
    method = "U-Pb"
    channels = Dict("d"=>"Pb207",
                    "D"=>"Pb206",
                    "P"=>"U238");
    standards = Dict("91500_zr" => "91500")
    glass = Dict("NIST612" => "NIST612")
    blk, fit = process!(myrun,method,channels,standards,glass,
                        nblank=2,ndrift=2,ndown=0)
    export2IsoplotR(myrun,method,channels,blk,fit;
                    fname="/home/pvermees/temp/NHM.json")
    results = averat(myrun,channels,blk,fit;method=method)
    CSV.write("/home/pvermees/temp/NHM.csv",results)
    p = KJ.plot(myrun[snum],method,channels,blk,fit,standards,glass;
                transformation="log",den=nothing)
elseif option == "NHM-carbonate"
    myrun = load("/home/pvermees/Documents/Plasmatrace/NHM/240408_ET_Carbonate.csv",
                 "/home/pvermees/Documents/Plasmatrace/NHM/240408_FP_CarbonateDating.Iolite.csv";
                 instrument="Agilent")
    method = "U-Pb"
    channels = Dict("d"=>"Pb207",
                    "D"=>"Pb206",
                    "P"=>"U238");
    standards = Dict("RA138_cc" => "RA138")
    glass = Dict("NIST612" => "NIST612")
    setBwin!(myrun[1],[(40,60)];seconds=true)
    blk, fit = process!(myrun,method,channels,standards,glass,
                        nblank=2,ndrift=2,ndown=0)
    export2IsoplotR(myrun,method,channels,blk,fit;
                    fname="/home/pvermees/temp/NHM-carbonate.json")
    results = averat(myrun,channels,blk,fit;method=method)
    CSV.write("/home/pvermees/temp/NHM-carbonate.csv",results)
    if true
        # rm("/home/pvermees/temp/NHM-carbonate.pdf")
        for samp in myrun
            savefig(KJ.plot(samp,method,channels,blk,fit,standards,glass;
                            transformation="log",den=nothing),
                    "/home/pvermees/temp/temp.pdf")
            append_pdf!("/home/pvermees/temp/NHM-carbonate.pdf",
                        "/home/pvermees/temp/temp.pdf",
                        cleanup=true)
        end
    end
elseif option == "KJgui"
    TUI(KJgui;reset=true)
else
    include("runtests.jl")
end
