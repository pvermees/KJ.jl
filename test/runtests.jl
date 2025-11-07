using KJ, Test, CSV, Infiltrator, DataFrames,
    Statistics, Distributions, Random
import Plots
include("helper.jl")

function loadtest(verbose=false)
    myrun = load("data/MWE";format="Agilent")
    if verbose summarise(myrun;verbose=true) end
    return myrun
end

function plottest(option="all")
    myrun = loadtest()
    if option in (1,"all")
        p = KJ.plot(myrun[1];
                    channels=["Hf176 -> 258","Hf178 -> 260"])
        @test display(p) != NaN
    end
    if option in (2,"all")
        p = KJ.plot(myrun[1];
                    channels=["Lu175 -> 175","Hf176 -> 258","Hf178 -> 260"],
                    den="Hf178 -> 260",
                    transformation = "log")
        @test display(p) != NaN
    end
end

function windowtest()
    myrun = loadtest()
    i = 2
    setSwin!(myrun[i],[(70,90),(100,140)])
    setBwin!(myrun[i],[(0,22)];seconds=true)
    setSwin!(myrun[i],[(37,65)];seconds=true)
    p = KJ.plot(myrun[i];channels=["Hf176 -> 258","Hf178 -> 260"])
    @test display(p) != NaN
end

function blanktest(;doplot=false,ylim=:auto,transformation=nothing)
    myrun = loadtest()
    blk = fitBlanks(myrun;nblank=2)
    if doplot
        p = KJ.plot(myrun[1],ylim=ylim,transformation=transformation)
        plotFittedBlank!(p,myrun[1],blk,transformation=transformation)
        @test display(p) != NaN
    end
    return myrun, blk
end

function mmediantest()
    n = 10
    v = collect(1:10)
    i = moving_median_indices(n;b=2)
    m = [median(v[i[j, :]]) for j in 1:n]
    println(v)
    display(i)
    println(m)
end

function outliertest()
    n = 100
    Random.seed!(0)
    random_values = rand(Distributions.Normal(0,1), n)
    random_values[50] = rand(Distributions.Normal(0,10),1)[1]
    outliers = detect_outliers(random_values)
    col = fill(1,n)
    col[outliers] .= 0
    p1 = Plots.scatter(1:n,random_values;
                       label=nothing,
                       marker_z=col,
                       legend=false)
    mu = [0.0; 2.0]
    Sigma = [1.0 -1.2; -1.2 2.0]
    mvn = Distributions.MvNormal(mu, Sigma)
    random_matrix = Matrix(rand(mvn, n)')
    #random_matrix[50,1] = rand(Distributions.Normal(3,10),1)[1]
    outliers_2 = detect_outliers(random_matrix)
    col = fill(1,n)
    col[outliers_2] .= 0
    p2 = Plots.scatter(random_matrix[:,1],
                       random_matrix[:,2];
                       label=nothing,
                       marker_z=col,
                       legend=false)
    p = Plots.plot(p1,p2;layout=(1,2))
    display(p)
end

function standardtest(verbose=false)
    myrun, blk = blanktest()
    standards = Dict("BP_gt" => "BP")
    setGroup!(myrun,standards)
    anchors = getStandardAnchors("Lu-Hf",standards)
    if verbose
        println(anchors)
        summarise(myrun;verbose=true,n=5)
    end
    return myrun
end

function chauvenetest()
    myrun = standardtest()
    channels = ["Hf178 -> 260","Hf176 -> 258","Lu175 -> 175"]
    chauvenet!(myrun;channels=channels)
end

function fixedLuHf(drift,down,mfrac,PAcutoff,adrift)
    myrun, blk = blanktest()
    method = "Lu-Hf"
    channels = Dict("d" => "Hf178 -> 260",
                    "D" => "Hf176 -> 258",
                    "P" => "Lu175 -> 175")
    glass = Dict("NIST612" => "NIST612p")
    setGroup!(myrun,glass)
    standards = Dict("BP_gt" => "BP")
    setGroup!(myrun,standards)
    fit = (drift=drift,down=down,mfrac=mfrac,
           PAcutoff=PAcutoff,adrift=adrift)
    return myrun, blk, method, channels, glass, standards, fit
end

function predictest()
    drift = [3.91]
    down = [0.0,0.0045]
    mfrac = -0.38
    myrun, blk, method, channels, glass, standards, fit =
        fixedLuHf(drift,down,mfrac,nothing,drift)
    samp = myrun[1]
    if samp.group == "sample"
        println("Not a standard")
        return samp, method, fit, blk, channels, standards, glass
    else
        pred = predict(samp,method,fit,blk,channels,standards,glass)
        p, offset = KJ.plot(samp,method,channels,blk,fit,standards,glass;
                            den="Hf176 -> 258",
                            transformation="log",
                            return_offset=true)
        @test display(p) != NaN
        return samp, method, fit, blk, channels, standards, glass, p, offset
    end
end
    
function partest(parname,paroffsetfact)
    samp,method,fit,blk,channels,standards,glass,p,offset = predictest()
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
        anchors = getStandardAnchors(method,standards)
        plotFitted!(p,samp,blk,adjusted_fit,channels,anchors;
                    transformation="log",offset=offset,
                    linecolor="red",debug=false)
    end
    @test display(p) != NaN
end

function driftest()
    partest("drift",1.0)
end

function downtest()
    partest("down",4.0)
end

function mfractest()
    partest("mfrac",0.2)
end

function fractionationtest(all=true)
    myrun, blk = blanktest()
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
        mf = fractionation(myrun,method,blk,channels,glass)
        fit = fractionation(myrun,method,blk,channels,standards,mf;
                            ndrift=1,ndown=1)
        println(fit)
        print("no glass: ")
        fit = fractionation(myrun,method,blk,channels,standards,nothing;
                            ndrift=1,ndown=1)
        println(fit)
        println("two joint steps: ")
    end
    fit = fractionation(myrun,"Lu-Hf",blk,channels,standards,glass;
                        ndrift=1,ndown=1)
    if (all)
        println(fit)
        return myrun, blk, fit, channels, standards, glass
    else
        Ganchors = getGlassAnchors(method,glass)
        Sanchors = getStandardAnchors(method,standards)
        anchors = merge(Sanchors,Ganchors)
        return myrun, blk, fit, channels, standards, glass, anchors
    end
end

function RbSrTest(show=true)
    myrun = load("data/Rb-Sr",format="Agilent")
    method = "Rb-Sr"
    channels = Dict("d"=>"Sr88 -> 104",
                    "D"=>"Sr87 -> 103",
                    "P"=>"Rb85 -> 85")
    standards = Dict("MDC_bt" => "MDC -")
    setGroup!(myrun,standards)
    blank = fitBlanks(myrun;nblank=2)
    fit = fractionation(myrun,method,blank,channels,standards,0.11935;
                        ndown=0,ndrift=1,verbose=false)
    anchors = getStandardAnchors(method,standards)
    if show
        p = KJ.plot(myrun[2],channels,blank,fit,anchors;
                    transformation="log",den="Sr87 -> 103")
        @test display(p) != NaN
    end
    export2IsoplotR(myrun,method,channels,blank,fit;
                    prefix="Entire",fname="output/Entire.json")
    return myrun, blank, fit, channels, standards, anchors
end

function KCaTest(show=true)
    myrun = load("data/K-Ca",format="Agilent")
    method = "K-Ca"
    channels = Dict("d"=>"Ca44 -> 63",
                    "D"=>"Ca40 -> 59",
                    "P"=>"K39 -> 39")
    standards = Dict("EntireCreek_bt" => "EntCrk")
    setGroup!(myrun,standards)
    blank = fitBlanks(myrun;nblank=2)
    fit = fractionation(myrun,method,blank,channels,standards,nothing;
                        ndown=0,ndrift=1,verbose=false)
    anchors = getStandardAnchors(method,standards)
    if show
        p = KJ.plot(myrun[3],channels,blank,fit,anchors,
                    transformation="log",den=nothing)
        @test display(p) != NaN
    end
    export2IsoplotR(myrun,method,channels,blank,fit;
                    prefix="EntCrk",fname="output/Entire_KCa.json")
    return myrun, blank, fit, channels, standards, anchors
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

function histest(;LuHf=false,show=true)
    if LuHf
        myrun,blk,fit,channels,standards,glass,anchors =
            fractionationtest(false)
        standard = "BP_gt"
    else
        myrun,blk,fit,channels,standards,anchors = RbSrTest(false)
        standard = "MDC_bt"
    end
    print(fit)
    dats, covs = pool(myrun;signal=true,group=standard,include_covmats=true)
    anchor = anchors[standard]
    pooled = DataFrame()
    pred = DataFrame()
    for i in eachindex(dats)
        pooled = vcat(pooled,dats[i])
        pred = vcat(pred,predict(dats[i],covs[i],fit,blk,channels,anchor))
    end
    Pm = pooled[:,channels["P"]]
    Dm = pooled[:,channels["D"]]
    dm = pooled[:,channels["d"]]
    Pp = pred[:,"P"]
    Dp = pred[:,"D"]
    dp = pred[:,"d"]
    if show
        plot_residuals(Pm,Dm,dm,Pp,Dp,dp)
        df = DataFrame(Pm=Pm,Dm=Dm,dm=dm,Pp=Pp,Dp=Dp,dp=dp)
        CSV.write("output/pooled_" * standard * ".csv",df)
    end
    return anchors, fit, Pm, Dm, dm
end

function processtest(show=true)
    myrun = load("data/MWE",format="Agilent")
    method = "Lu-Hf";
    channels = Dict("d"=>"Hf178 -> 260",
                    "D"=>"Hf176 -> 258",
                    "P"=>"Lu175 -> 175")
    standards = Dict("Hogsbo_gt" => "hogsbo")
    glass = Dict("NIST612" => "NIST612p")
    blk, fit = process!(myrun,method,channels,standards,glass;
                        nblank=2,ndrift=1,ndown=1,verbose=false)
    if show
        p = KJ.plot(myrun[2],method,channels,blk,fit,standards,glass;
                    transformation="log",den="Hf176 -> 258")
        @test display(p) != NaN
    end
    return myrun, method, channels, blk, fit
end

function PAtest(verbose=false)
    myrun = load("data/Lu-Hf",format="Agilent")
    method = "Lu-Hf"
    channels = Dict("d"=>"Hf178 -> 260",
                    "D"=>"Hf176 -> 258",
                    "P"=>"Lu175 -> 175")
    standards = Dict("Hogsbo_gt" => "hogsbo")
    glass = Dict("NIST612" => "NIST612p")
    cutoff = 1e7
    blk, fit = process!(myrun,method,channels,standards,glass;
                        PAcutoff=cutoff,nblank=2,ndrift=1,ndown=1)
    ratios = averat(myrun,channels,blk,fit)
    if verbose println(first(ratios,5)) end
    return ratios
end

function exporttest()
    ratios = PAtest()
    selection = prefix2subset(ratios,"BP") # "hogsbo"
    CSV.write("output/BP.csv",selection)
    export2IsoplotR(selection,"Lu-Hf",fname="output/BP.json")
end

function UPbtest()
    myrun = load("data/U-Pb",format="Agilent",head2name=false)
    method = "U-Pb"
    standards = Dict("Plesovice_zr" => "STDCZ",
                     "91500_zr" => "91500")
    glass = Dict("NIST610" => "610",
                 "NIST612" => "612")
    channels = Dict("d"=>"Pb207","D"=>"Pb206","P"=>"U238")
    blank, fit = process!(myrun,"U-Pb",channels,standards,glass;
                          nblank=2,ndrift=1,ndown=1)
    export2IsoplotR(myrun,method,channels,blank,fit;
                    fname="output/UPb.json")
    p = KJ.plot(myrun[37],method,channels,blank,fit,standards,glass;
                transformation="log",den="Pb206")
    @test display(p) != NaN
end

function iCaptest(verbose=true)
    myrun = load("data/iCap",format="ThermoFisher")
    if verbose summarise(myrun;verbose=true,n=5) end
end

function carbonatetest(verbose=false)
    method = "U-Pb"
    myrun = load("data/carbonate",format="Agilent")
    standards = Dict("WC1_cc"=>"WC1")
    glass = Dict("NIST612"=>"NIST612")
    channels = Dict("d"=>"Pb207","D"=>"Pb206","P"=>"U238")
    blk, fit = process!(myrun,method,channels,standards,glass;
                        nblank=2,ndrift=1,ndown=1,verbose=verbose)
    export2IsoplotR(myrun,method,channels,blk,fit;
                    prefix="Duff",fname="output/Duff.json")
    p = KJ.plot(myrun[4],method,channels,blk,fit,standards,glass;
                transformation="log")
    @test display(p) != NaN
end

function timestamptest(verbose=true)
    myrun = load("data/timestamp/Moreira_data.csv",
                 "data/timestamp/Moreira_timestamps.csv";
                 format="Agilent")
    if verbose summarise(myrun;verbose=true,n=5) end
    p = KJ.plot(myrun[2];
                transformation="sqrt")
    @test display(p) != NaN
end

function mineraltest()
    internal = getInternal("zircon","Si29")
end

function concentrationtest()
    method = "concentrations"
    myrun = load("data/Lu-Hf",format="Agilent")
    internal = ("Al27 -> 27",1.2e5)
    glass = Dict("NIST612" => "NIST612p")
    setGroup!(myrun,glass)
    blk, fit = process!(myrun,internal,glass;nblank=2)
    p = KJ.plot(myrun[4],blk,fit,internal[1];
                transformation="log",den=internal[1])
    conc = concentrations(myrun,blk,fit,internal)
    @test display(p) != NaN
    return conc
end

function internochrontest(show=true)
    myrun = load("data/lines",format="Agilent")
    method = "Lu-Hf"
    channels = Dict("d"=>"Hf178 -> 260",
                    "D"=>"Hf176 -> 258",
                    "P"=>"Lu175 -> 175")
    standards = Dict("Hogsbo_gt" => "Hog",
                     "BP_gt" => "BP")
    glass = Dict("NIST610" => "NIST610")
    blk, fit = process!(myrun,method,channels,standards,glass)
    isochron = internochron(myrun,channels,blk,fit;method=method)
    CSV.write("output/isochron.csv",isochron)
    if show
        p = internoplot(myrun[11],channels,blk,fit;method=method)
        @test display(p) != NaN
    end
end

function internochronUPbtest(show=true)
    method = "U-Pb"
    myrun = load("data/carbonate",format="Agilent")
    standards = Dict("WC1_cc"=>"WC1")
    glass = Dict("NIST612"=>"NIST612")
    channels = Dict("d"=>"Pb207","D"=>"Pb206","P"=>"U238")
    blk, fit = process!(myrun,method,channels,standards,glass,
                        nblank=2,ndrift=1,ndown=1)
    isochron = internochron(myrun,channels,blk,fit;method=method)
    CSV.write("output/isochronUPb.csv",isochron)
    if show
        p = internoplot(myrun[7],channels,blk,fit;method=method)
        @test display(p) != NaN
    end
end

function maptest()
    method = "concentrations"
    myrun = load("data/timestamp/NHM_cropped.csv",
                 "data/timestamp/NHM_timestamps.csv";
                 format="Agilent")
    internal = getInternal("zircon","Si29")
    glass = Dict("NIST612" => "NIST612")
    setGroup!(myrun,glass)
    blk, fit = process!(myrun,internal,glass;nblank=2)
    conc = concentrations(myrun[10],blk,fit,internal)
    p = plotMap(conc,"ppm[U] from U238";
                clims=(0,500),
                colorbar_scale=:identity)
    CSV.write("output/Umap.csv",conc)
    @test display(p) != NaN
end

function map_dating_test()
    method = "U-Pb"
    myrun = load("data/timestamp/NHM_cropped.csv",
                 "data/timestamp/NHM_timestamps.csv";
                 format="Agilent")
    standards = Dict("91500_zr"=>"91500")
    glass = Dict("NIST612" => "NIST612")
    channels = Dict("d"=>"Pb207","D"=>"Pb206","P"=>"U238")
    blk, fit = process!(myrun,method,channels,standards,glass,
                        nblank=2,ndrift=1,ndown=0)
    snum = 10
    P,D,d,x,y = atomic(myrun[snum],channels,blk,fit;add_xy=true)
    df = DataFrame(P=P,D=D,d=d,x=x,y=y)
    p = plotMap(df,"P";clims=(1e3,1e6))
    @test display(p) != NaN
end

function map_fail_test()
    myrun, method, channels, blk, fit = processtest()
    P,D,d,x,y = atomic(myrun[2],channels,blk,fit;add_xy=true)
    df = DataFrame(P=P,D=D,d=d,x=x,y=y)
    p = plotMap(df,"P")
    @test display(p) != NaN
    conc = concentrationtest()
    p = plotMap(conc,"Lu175 -> Lu175")
    @test display(p) != NaN
end

function glass_only_test()
    myrun = load("data/U-Pb",format="Agilent",head2name=false)
    method = "U-Pb"
    glass = standards = Dict("NIST610" => "610",
                             "NIST612" => "612")
    channels = Dict("d"=>"Pb207","D"=>"Pb206","P"=>"U238")
    blank, pars = process!(myrun,"U-Pb",channels,standards,glass;
                           nblank=2,ndrift=1,ndown=1)
    export2IsoplotR(myrun,method,channels,blank,pars;
                    fname="output/UPb_with_glass.json")
    p = KJ.plot(myrun[37],method,channels,blank,pars,standards,glass;
                transformation="log",den="Pb206")
    display(p)
end

function synthetictest(truefit;kw...)
    method = "Lu-Hf"
    channels = Dict("d"=>"Hf178 -> 260",
                    "D"=>"Hf176 -> 258",
                    "P"=>"Lu175 -> 175")
    standards = Dict("BP_gt" => "BP")
    glass = Dict("NIST612" => "NIST612")
    myrun, channels = synthetic(;
                                truefit=truefit,
                                lambda=1.867e-5,
                                t_std=1745.0,
                                y0_std=3.55,
                                t_smp=1029.7,
                                y0_smp=3.55,
                                y0_glass=3.544842,
                                channels=channels,
                                kw...)
    return method, channels, standards, glass, myrun
end

function SStest()
    truefit = (drift=[0.0],down=[0.0],mfrac=0.0,
               PAcutoff=nothing,adrift=[0.0])
    x0_std = 1.0
    y0_std = 1.0
    D_blank = 1.0
    down = truefit.down
    drift = truefit.drift
    mfrac = truefit.mfrac
    method, channels, standards, glass, myrun = synthetictest(truefit)
    blk, fit = process!(myrun,method,channels,standards,glass;
                        nblank=2,ndrift=length(drift),
                        ndown=length(down),verbose=false)
    nstep = 20
    ss = fill(0.0,nstep)
    fit = deepcopy(truefit)
    dd = range(start=fit.drift[1]-1.0,stop=fit.drift[1]+1.0,length=nstep)
    for i in eachindex(dd)
        fit.drift[1] = dd[i]
        ss[i] = SS(fit,myrun,method,standards,blk,channels)
    end
    p = Plots.plot(dd,ss,seriestype=:line,label="drift")
    dwn = range(start=fit.down[1]-1.0,stop=fit.down[1]+1.0,length=nstep)
    fit = deepcopy(truefit)
    for i in eachindex(dwn)
        fit.down[1] = dwn[i]
        ss[i] = SS(fit,myrun,method,standards,blk,channels)
    end
    Plots.plot!(p,dwn,ss,seriestype=:line,linecolor=:red,label="down")
    SS_truth = SS(truefit,myrun,method,standards,blk,channels;verbose=true)
    display(p)
end

function accuracytest(;drift=[0.0],down=[0.0],mfrac=0.0,show=true,kw...)
    truefit=(drift=drift,down=down,mfrac=mfrac,PAcutoff=nothing,adrift=drift)
    method, channels, standards, glass, myrun = synthetictest(truefit,kw...)
    blk, fit = process!(myrun,method,channels,standards,glass;
                        nblank=2,ndrift=1,ndown=1,verbose=false)
    SS_solution = SS(fit,myrun,method,standards,blk,channels)
    SS_truth = SS(truefit,myrun,method,standards,blk,channels)
    @test SS_solution < SS_truth
    if show
        den = nothing # "Hf176 -> 258" #
        p1 = KJ.plot(myrun[1],method,channels,blk,fit,standards,glass;
                     transformation="sqrt",den=den)
        p2 = KJ.plot(myrun[3],method,channels,blk,fit,standards,glass;
                     transformation="sqrt",den=den)
        p3 = KJ.plot(myrun[4],method,channels,blk,fit,standards,glass;
                     transformation="sqrt",den=den)
        p = Plots.plot(p1,p2,p3,layout=(1,3))
        @test display(p) != NaN
    end
end

module test
function extend!(_KJ::AbstractDict)
    old = _KJ["tree"]["top"]
    _KJ["tree"]["top"] = (message = "test",
                          help = "test",
                          action = old.action)
end
export KJtree!
end
using .test
function extensiontest(verbose=true)
    TUI(test;logbook="logs/extension.log")
end

function TUItest()
    TUI(;logbook="logs/Lu-Hf.log",reset=true)
end

Plots.closeall()

if true
    #=@testset "load" begin loadtest(true) end
    @testset "plot raw data" begin plottest(2) end
    @testset "set selection window" begin windowtest() end
    @testset "set method and blanks" begin blanktest() end
    @testset "moving median test" begin mmediantest() end=#
    @testset "outlier detection" begin outliertest() end
    #=@testset "assign standards" begin standardtest(true) end
    @testset "chauvenet test" begin chauvenetest() end
    @testset "predict" begin predictest() end
    @testset "predict drift" begin driftest() end
    @testset "predict down" begin downtest() end
    @testset "predict mfrac" begin mfractest() end
    @testset "fractionation" begin fractionationtest(true) end
    @testset "Rb-Sr" begin RbSrTest() end
    @testset "K-Ca" begin KCaTest() end
    @testset "hist" begin histest() end
    @testset "process run" begin processtest() end
    @testset "PA test" begin PAtest(true) end
    @testset "export" begin exporttest() end
    @testset "U-Pb" begin UPbtest() end
    @testset "iCap" begin iCaptest() end
    @testset "carbonate" begin carbonatetest() end
    @testset "timestamp" begin timestamptest() end
    @testset "stoichiometry" begin mineraltest() end
    @testset "concentration" begin concentrationtest() end
    @testset "Lu-Hf internochron" begin internochrontest() end
    @testset "UPb internochron" begin internochronUPbtest() end
    @testset "concentration map" begin maptest() end
    @testset "isotope ratio map" begin map_dating_test() end
    @testset "map fail test" begin map_fail_test() end
    @testset "glass as age standard test" begin glass_only_test() end
    @testset "extension test" begin extensiontest() end
    @testset "synthetic data" begin SStest() end
    @testset "accuracy test 1" begin accuracytest() end
    @testset "accuracy test 2" begin accuracytest(drift=[-2.0]) end
    @testset "accuracy test 3" begin accuracytest(mfrac=2.0) end
    @testset "accuracy test 4" begin accuracytest(down=[0.0,0.5]) end
    @testset "TUI test" begin TUItest() end=#

else
    TUI()
end
