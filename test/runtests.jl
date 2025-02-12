using KJ, Test, CSV, Infiltrator, DataFrames, Statistics
import Plots

function loadtest(verbose=false)
    myrun = load("data/Lu-Hf";instrument="Agilent")
    if verbose summarise(myrun;verbose=true,n=5) end
    return myrun
end

function plottest()
    myrun = loadtest()
    p = plot(myrun[1],["Hf176 -> 258","Hf178 -> 260"])
    @test display(p) != NaN
    p = plot(myrun[1],["Hf176 -> 258","Hf178 -> 260"], den="Hf178 -> 260")
    @test display(p) != NaN
end

function windowtest()
    myrun = loadtest()
    i = 2
    setSwin!(myrun[i],[(70,90),(100,140)])
    setBwin!(myrun[i],[(0,22)];seconds=true)
    setSwin!(myrun[i],[(37,65)];seconds=true)
    p = plot(myrun[i],["Hf176 -> 258","Hf178 -> 260"])
    @test display(p) != NaN
end

function blanktest(poisson=true)
    myrun = loadtest()
    blk = fitBlanks(myrun;nblank=2)
    dt = poisson ? dwelltime(myrun) : nothing
    return myrun, dt, blk
end

function standardtest(verbose=false)
    myrun, dt, blk = blanktest()
    standards = Dict("BP_gt" => "BP")
    setGroup!(myrun,standards)
    anchors = getAnchors("Lu-Hf",standards)
    if verbose
        println(anchors)
        summarise(myrun;verbose=true,n=5)
    end
end

function fixedLuHf(poisson=true)
    myrun, dt, blk = blanktest(poisson)
    method = "Lu-Hf"
    channels = Dict("d" => "Hf178 -> 260",
                    "D" => "Hf176 -> 258",
                    "P" => "Lu175 -> 175")
    glass = Dict("NIST612" => "NIST612p")
    setGroup!(myrun,glass)
    standards = Dict("BP_gt" => "BP")
    setGroup!(myrun,standards)
    fit = (drift=[-3.9225],
           down=[0.0,0.03362],
           mfrac=0.38426,
           PAcutoff=nothing,
           adrift=[-3.9225])
    return myrun, dt, blk, method, channels, glass, standards, fit
end

function predictest(poisson=true)
    myrun, dt, blk, method, channels, glass, standards, fit = fixedLuHf(poisson)
    samp = myrun[105]
    if samp.group == "sample"
        println("Not a standard")
    else
        pred = predict(samp,method,fit,blk,channels,standards,glass;dt=dt)
        p = plot(samp,method,channels,blk,fit,standards,glass;
                 dt=dt,transformation="log")
        @test display(p) != NaN
    end
    return samp,dt,method,fit,blk,channels,standards,glass,p
end
    
function partest(parname,paroffsetfact;poisson=true)
    samp,dt,method,fit,blk,channels,standards,glass,p = predictest(poisson)
    drift = fit.drift[1]
    down = fit.down[2]
    mfrac = fit.mfrac[1]
    for paroffset in paroffsetfact .* [-1,1]
        if parname=="drift"
            drift = fit.drift[1] + paroffset
        elseif parname=="down"
            down = fit.down[2] + paroffset
        elseif parname=="mfrac"
            mfrac = fit.mfrac[1] + paroffset
        end
        adjusted_fit = (drift=[drift],
                        down=[0.0,down],
                        mfrac=mfrac,
                        PAcutoff=nothing,
                        adrift=[drift])
        anchors = getAnchors(method,standards,false)
        offset = getOffset(samp,channels,blk,adjusted_fit,anchors,"log";
                           dt=dt)
        plotFitted!(p,samp,blk,adjusted_fit,channels,anchors;
                    offset=offset,transformation="log",linecolor="red",dt=dt)
    end
    @test display(p) != NaN
end

function driftest(poisson=true)
    partest("drift",1.0;poisson=poisson)
end

function downtest(poisson=true)
    partest("down",4.0;poisson=poisson)
end

function mfractest(poisson=true)
    partest("mfrac",0.2;poisson=poisson)
end

function fractionationtest(all=true;poisson=true)
    myrun, dt, blk = blanktest(poisson)
    method = "Lu-Hf"
    channels = Dict("d" => "Hf178 -> 260",
                    "D" => "Hf176 -> 258",
                    "P" => "Lu175 -> 175")
    glass = Dict("NIST612" => "NIST612p")
    setGroup!(myrun,glass)
    standards = Dict("BP_gt" => "BP")
    setGroup!(myrun,standards)
    if all
        println("two separate steps: ")
        mf = fractionation(myrun,method,blk,channels,glass;dt=dt)
        fit = fractionation(myrun,method,blk,channels,standards,mf;
                            dt=dt,ndrift=1,ndown=1)
        println(fit)
        print("no glass: ")
        fit = fractionation(myrun,method,blk,channels,standards,nothing;
                            ndrift=1,ndown=1,dt=dt)
        println(fit)
        println("two joint steps: ")
    end
    fit = fractionation(myrun,"Lu-Hf",blk,channels,standards,glass;
                        ndrift=1,ndown=1,dt=dt)
    if (all)
        println(fit)
        return myrun, dt, blk, fit, channels, standards, glass
    else
        Ganchors = getAnchors(method,glass,true)
        Sanchors = getAnchors(method,standards,false)
        anchors = merge(Sanchors,Ganchors)
        return myrun, dt, blk, fit, channels, standards, glass, anchors
    end
end

function RbSrTest(show=true;poisson=true)
    myrun = load("data/Rb-Sr",instrument="Agilent")
    method = "Rb-Sr"
    channels = Dict("d"=>"Sr88 -> 104",
                    "D"=>"Sr87 -> 103",
                    "P"=>"Rb85 -> 85")
    standards = Dict("MDC_bt" => "MDC -")
    setGroup!(myrun,standards)
    blank = fitBlanks(myrun;nblank=2)
    dt = dwelltime(myrun;poisson=poisson)
    fit = fractionation(myrun,method,blank,channels,standards,8.37861;
                        ndown=0,ndrift=1,dt=dt,verbose=false)
    anchors = getAnchors(method,standards)
    if show
        p = plot(myrun[2],channels,blank,fit,anchors,
                 transformation="log")#,den="Sr88 -> 104";dt=dt)
        @test display(p) != NaN
    end
    export2IsoplotR(myrun,method,channels,blank,fit;
                    prefix="Entire",fname="Entire.json",dt=dt)
    return myrun, dt, blank, fit, channels, standards, anchors
end

function KCaTest(show=true;poisson=true)
    myrun = load("data/K-Ca",instrument="Agilent")
    method = "K-Ca"
    channels = Dict("d"=>"Ca44 -> 63",
                    "D"=>"Ca40 -> 59",
                    "P"=>"K39 -> 39")
    standards = Dict("EntireCreek_bt" => "EntCrk")
    setGroup!(myrun,standards)
    blank = fitBlanks(myrun;nblank=2)
    dt = dwelltime(myrun;poisson=poisson)
    fit = fractionation(myrun,method,blank,channels,standards,nothing;
                        ndown=0,ndrift=1,dt=dt,verbose=false)
    anchors = getAnchors(method,standards)
    if show
        p = plot(myrun[3],channels,blank,fit,anchors,
                 transformation="log",den=nothing;dt=dt)
        @test display(p) != NaN
    end
    export2IsoplotR(myrun,method,channels,blank,fit,
                    prefix="EntCrk",fname="Entire_KCa.json";dt=dt)
    return myrun, dt, blank, fit, channels, standards, anchors
end

function plot_residuals(Pm,Dm,dm,Pp,Dp,dp)
    Pmisfit = Pm.-Pp
    Dmisfit = Dm.-Dp
    dmisfit = dm.-dp
    pP = Plots.histogram(Pmisfit,xlab="εP")
    pD = Plots.histogram(Dmisfit,xlab="εD")
    pd = Plots.histogram(dmisfit,xlab="εd")
    pPD = Plots.plot(Pmisfit,Dmisfit;
                     seriestype=:scatter,xlab="εP",ylab="εD")
    pPd = Plots.plot(Pmisfit,dmisfit;
                     seriestype=:scatter,xlab="εP",ylab="εd")
    pDd = Plots.plot(Dmisfit,dmisfit;
                     seriestype=:scatter,xlab="εD",ylab="εd")
    pDP = Plots.plot(Dmisfit,Pmisfit;
                     seriestype=:scatter,xlab="εD",ylab="εP")
    pdP = Plots.plot(dmisfit,Pmisfit;
                     seriestype=:scatter,xlab="εd",ylab="εP")
    pdD = Plots.plot(dmisfit,Dmisfit;
                     seriestype=:scatter,xlab="εd",ylab="εD")
    p = Plots.plot(pP,pDP,pdP,
                   pPD,pD,pdD,
                   pPd,pDd,pd;legend=false)
    @test display(p) != NaN    
end

function histest(;LuHf=false,show=true,poisson=true)
    if LuHf
        myrun,dt,blk,fit,channels,standards,glass,anchors =
            fractionationtest(false;poisson=poisson)
        standard = "BP_gt"
    else
        myrun,dt,blk,fit,channels,standards,anchors = RbSrTest(false;poisson=poisson)
        standard = "MDC_bt"
    end
    print(fit)
    pooled = pool(myrun;signal=true,group=standard)
    anchor = anchors[standard]
    pred = predict(pooled,fit,blk,channels,anchor;dt=dt)
    Pm = pooled[:,channels["P"]]
    Dm = pooled[:,channels["D"]]
    dm = pooled[:,channels["d"]]
    Pp = pred[:,"P"]
    Dp = pred[:,"D"]
    dp = pred[:,"d"]
    if show
        plot_residuals(Pm,Dm,dm,Pp,Dp,dp)
        df = DataFrame(Pm=Pm,Dm=Dm,dm=dm,Pp=Pp,Dp=Dp,dp=dp)
        CSV.write("pooled_" * standard * ".csv",df)
    end
    return anchors, fit, Pm, Dm, dm
end

function processtest(poisson=true)
    myrun = load("data/Lu-Hf",instrument="Agilent")
    dt = poisson ? dwelltime(myrun) : nothing
    method = "Lu-Hf";
    channels = Dict("d"=>"Hf178 -> 260",
                    "D"=>"Hf176 -> 258",
                    "P"=>"Lu175 -> 175");
    standards = Dict("Hogsbo_gt" => "hogsbo")#"BP_gt" => "BP")#
    glass = Dict("NIST612" => "NIST612p")
    blk, fit = process!(myrun,method,channels,standards,glass;
                        nblank=2,ndrift=1,ndown=1,dt=dt)
    p = plot(myrun[2],method,channels,blk,fit,standards,glass;
             transformation="log",den=nothing,dt=dt)
    @test display(p) != NaN
end

function PAtest(verbose=false;poisson=true)
    myrun = load("data/Lu-Hf",instrument="Agilent")
    dt = poisson ? dwelltime(myrun) : nothing
    method = "Lu-Hf"
    channels = Dict("d"=>"Hf178 -> 260",
                    "D"=>"Hf176 -> 258",
                    "P"=>"Lu175 -> 175")
    standards = Dict("Hogsbo_gt" => "hogsbo")
    glass = Dict("NIST612" => "NIST612p")
    cutoff = 1e7
    blk, fit = process!(myrun,method,channels,standards,glass;
                        PAcutoff=cutoff,nblank=2,ndrift=1,ndown=1,dt=dt)
    ratios = averat(myrun,channels,blk,fit;method=method,dt=dt)
    if verbose println(first(ratios,5)) end
    return ratios
end

function exporttest(poisson=true)
    ratios = PAtest(poisson=poisson)
    prefix = "BP" # hogsbo
    selection = prefix2subset(ratios,prefix)
    fname = poisson ? "poisson" * prefix : "non-poisson" * prefix
    CSV.write(fname * ".csv",selection)
    export2IsoplotR(selection,"Lu-Hf",fname=prefix * ".json")
end

function UPbtest(poisson=true)
    myrun = load("data/U-Pb",instrument="Agilent",head2name=false)
    dt = poisson ? dwelltime(myrun) : nothing
    method = "U-Pb"
    standards = Dict("Plesovice_zr" => "STDCZ",
                     "91500_zr" => "91500")
    glass = Dict("NIST610" => "610",
                 "NIST612" => "612")
    channels = Dict("d"=>"Pb207","D"=>"Pb206","P"=>"U238")
    blank, pars = process!(myrun,"U-Pb",channels,standards,glass,
                           nblank=2,ndrift=1,ndown=1,dt=dt)
    export2IsoplotR(myrun,method,channels,blank,pars;
                    fname="UPb.json",dt=dt)
    p = plot(myrun[37],method,channels,blank,pars,standards,glass;
             transformation="log",den="Pb206",dt=dt)
    @test display(p) != NaN
end

function iCaptest(verbose=true)
    myrun = load("data/iCap",instrument="ThermoFisher")
    if verbose summarise(myrun;verbose=true,n=5) end
end

function carbonatetest(verbose=false;poisson=true)
    method = "U-Pb"
    myrun = load("data/carbonate",instrument="Agilent")
    dt = poisson ? dwelltime(myrun) : nothing
    standards = Dict("WC1_cc"=>"WC1")
    glass = Dict("NIST612"=>"NIST612")
    channels = Dict("d"=>"Pb207","D"=>"Pb206","P"=>"U238")
    blk, fit = process!(myrun,method,channels,standards,glass;
                        nblank=2,ndrift=1,ndown=1,verbose=verbose,dt=dt)
    export2IsoplotR(myrun,dt,method,channels,blk,fit;
                    prefix="Duff",fname="Duff.json")
    p = plot(myrun[3],method,channels,blk,fit,standards,glass;
             transformation=nothing,num=["Pb207"],
             den="Pb206",ylim=[-0.02,0.3],dt=dt)
    @test display(p) != NaN
end

function timestamptest(verbose=true)
    myrun = load("data/timestamp/MSdata.csv",
                 "data/timestamp/timestamp.csv";
                 instrument="Agilent")
    if verbose summarise(myrun;verbose=true,n=5) end
    p = plot(myrun[2];transformation="sqrt")
    @test display(p) != NaN
end

function mineraltest()
    internal = getInternal("zircon","Si29")
end

function concentrationtest()
    method = "concentrations"
    myrun = load("data/Lu-Hf",instrument="Agilent")
    internal = ("Al27 -> 27",1.2e5)
    glass = Dict("NIST612" => "NIST612p")
    setGroup!(myrun,glass)
    blk, fit = process!(myrun,internal,glass;nblank=2)
    p = plot(myrun[4],blk,fit,internal[1];
             transformation="log",den=internal[1])
    conc = concentrations(myrun,blk,fit,internal)
    @test display(p) != NaN
end

module test
function extend!(_KJ::AbstractDict)
    old = _KJ["tree"]["top"]
    _KJ["tree"]["top"] = (message = "test", help = "test", action = old.action)
end
export KJtree!
end
using .test
function extensiontest(verbose=true)
    TUI(test,logbook="logs/extension.log")
end

function TUItest()
    TUI(logbook="logs/test.log",reset=true)
end

function GUItest()
    TUI(KJgui;logbook="logs/KJgui.log",reset=true)
end

Plots.closeall()

if true
    @testset "load" begin loadtest(true) end
    @testset "plot raw data" begin plottest() end
    @testset "set selection window" begin windowtest() end
    @testset "poisson blank test" begin blanktest() end
    @testset "non-poisson blank test" begin blanktest(false) end
    @testset "assign standards" begin standardtest(true) end
    @testset "poisson predict" begin predictest() end
    @testset "non-poisson predict" begin predictest(false) end
    @testset "poisson drift test" begin driftest() end
    @testset "non-poisson drift test" begin driftest(false) end
    @testset "poisson downhole test" begin downtest() end
    @testset "non-poisson downhole test" begin downtest(false) end
    @testset "poisson mass fractionation" begin mfractest() end
    @testset "non-poisson mass fractionation" begin mfractest(false) end
    @testset "poisson fractionation" begin fractionationtest(true) end
    @testset "non-poisson fractionation" begin fractionationtest(true;poisson=false) end
    @testset "poisson Rb-Sr" begin RbSrTest() end
    @testset "non-poisson Rb-Sr" begin RbSrTest(poisson=false) end
    @testset "poisson K-Ca" begin KCaTest() end
    @testset "non-poisson K-Ca" begin KCaTest(poisson=false) end
    @testset "poisson hist" begin histest() end
    @testset "non-poisson hist" begin histest(poisson=false) end
    @testset "poisson process" begin processtest() end
    @testset "non-poisson process" begin processtest(false) end
    @testset "poisson PA" begin PAtest(true) end
    @testset "non-poisson PA" begin PAtest(true;poisson=false) end
    @testset "poisson export" begin exporttest(false) end
    @testset "non-poisson export" begin exporttest(true) end
    @testset "poisson U-Pb" begin UPbtest() end
    @testset "non-poisson U-Pb" begin UPbtest(false) end
    @testset "iCap test" begin iCaptest() end
    @testset "poisson carbonate" begin carbonatetest() end
    @testset "non-poisson carbonate" begin carbonatetest(poisson=false) end
    @testset "timestamp" begin timestamptest() end
    @testset "stoichiometry" begin mineraltest() end
    @testset "concentration" begin concentrationtest() end
    @testset "extension" begin extensiontest() end
    @testset "TUI" begin TUItest() end
    # @testset "KJgui" begin GUItest() end
else
    TUI()
end
