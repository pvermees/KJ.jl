using Test, KJ, Dates, DataFrames, Infiltrator
import Plots, Distributions, CSV, Statistics,
    Random, LinearAlgebra, Aqua
include("synthetic.jl")

function loadtest(;dname="data/Lu-Hf",
                  verbose=false)
    myrun = load(dname;format="Agilent")
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

function windowtest(show=true)
    myrun = load("data/Lu-Hf";format="Agilent")
    setBwin!(myrun)
    setSwin!(myrun)
    i = 2
    setSwin!(myrun[i],[(70,90),(100,140)])
    setBwin!(myrun[i],[(0,22)];seconds=true)
    setSwin!(myrun[i],[(37,65)];seconds=true)
    if show
        p = KJ.plot(myrun[i];channels=["Hf176 -> 258","Hf178 -> 260"])
        @test display(p) != NaN
    end
    return myrun
end

function blanktest(;myrun=loadtest(),
                   doplot=false,
                   ylim=:auto,
                   transformation=nothing)
    blk = fitBlanks(myrun;nblank=2)
    if doplot
        p = KJ.plot(myrun[1];ylim=ylim,transformation=transformation)
        plotFittedBlank!(p,myrun[1],blk,transformation=transformation)
        @test display(p) != NaN
    end
    return blk
end

function mmediantest()
    n = 10
    v = collect(1:10)
    i = moving_median_indices(n;b=2)
    m = [Distributions.median(v[i[j, :]]) for j in 1:n]
    println(v)
    display(i)
    println(m)
end

function outliertest_synthetic()
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
    random_matrix[50,1] = rand(Distributions.Normal(3,10),1)[1]
    outliers_2 = detect_outliers(random_matrix)
    col = fill(1,n)
    col[outliers_2] .= 0
    p2 = Plots.scatter(random_matrix[:,1],
                       random_matrix[:,2];
                       label=nothing,
                       marker_z=col,
                       legend=false)
    p = Plots.plot(p1,p2;layout=(1,2))
    @test display(p) != NaN
end

function outliertest_sample(show=true)
    myrun = load("data/K-Ca";format="Agilent")
    setBwin!(myrun)
    setSwin!(myrun)
    channels = ["K39 -> 39","Ca40 -> 59","Ca44 -> 63"]
    detect_outliers!(myrun;include_samples=true,channels=channels)
    if show
        p = KJ.plot(myrun[1];
                    channels=channels,
                    transformation="log",
                    den="Ca40 -> 59")
        @test display(p) != NaN
    end
end

function methodtest()
    interferences = Interference(ions = Dict("Hf176" => ["Lu176","Yb176"]),
                                 proxies = Dict("Lu176" => "Lu175","Yb176" => "Yb172"),
                                 channels = Dict("Lu175" => "Lu175 -> 257","Yb172" => "Yb172 -> 254"),
                                 bias = Dict("Lu175" => ["some_Lu_standard"], "Yb172" => ["some_REE_standard"]))
    method = Gmethod("Lu-Hf";
                     channels=(P="Lu175 -> 175",D="Hf176 -> 258",d="Hf178 -> 260"),
                     standards=["Hogsbo"],
                     bias=Dict("Hf" => ["NIST612"]),
                     interference=interferences)
    @test method.fractionation.proxies.d == "Hf178"
end

function refmattest(verbose=false)
    myrun = loadtest()
    refmats = Dict("BP" => "BP")
    setGroup!(myrun,refmats)
    if verbose
        summarise(myrun;verbose=true,n=5)
    end
    return myrun
end

function getmethod(name::AbstractString,
                   refmats::AbstractDict)
    standards = collect(keys(refmats))
    if name == "U-Pb"
        return Gmethod(name;standards=standards)
    end
    ndrift = 2
    ndown = 1
    p = nothing
    c = nothing
    if name=="Lu-Hf"
        p = (P="Lu175",D="Hf176",d="Hf178")
        c = (P="Lu175 -> 175",D="Hf176 -> 258",d="Hf178 -> 260")
    elseif name=="Rb-Sr"
        p = (P="Rb85",D="Sr87",d="Sr88")
        c = (P="Rb85 -> 85",D="Sr87 -> 103",d="Sr88 -> 104")
    elseif name=="K-Ca"
        p = (P="K39",D="Ca40",d="Ca44")
        c = (P="K39 -> 39",D="Ca40 -> 59",d="Ca44 -> 63")
        ndrift = 1
        ndown = 0
    end
    return Gmethod(name;
                   standards=standards,
                   proxies=p,channels=c,
                   ndrift=ndrift,ndown=ndown)
end

function predictsettings(option::AbstractString="Lu-Hf")
    dname = nothing
    refmats = nothing
    drift = nothing
    down = nothing
    if option=="Lu-Hf"
        dname = "data/Lu-Hf"
        refmats = Dict("BP" => "BP")
        drift = [4.22,0.0]
        down = [0.0,0.15]
    elseif option=="Rb-Sr"
        dname = "data/Rb-Sr"
        refmats = Dict("MDC" => "MDC -")
        drift = [1.0,0.0]
        down = [0.0,0.14]
    elseif option=="K-Ca"
        dname = "data/K-Ca"
        refmats = Dict("MDC" => "MDC_")
        drift = [100.0]
        down = [0.0]
    end
    myrun = loadtest(;dname=dname)
    method = getmethod(option,refmats)
    setGroup!(myrun,refmats)
    detect_outliers!(myrun;channels=getChannels(method))
    method.PAcutoff = nothing
    fit = Gfit(method;drift=drift,down=down)
    return myrun, method, fit
end

function predictest(option="Lu-Hf";
                    snum=1,
                    transformation="log")
    myrun, method, fit = predictsettings(option)
    fit.blank = blanktest(;myrun=myrun)
    samp = myrun[snum]
    if samp.group == "sample"
        println("Not a standard")
        return samp, method, fit
    else
        p, offset = KJ.plot(samp,method;fit=fit,
                            den=method.fractionation.channels.D,
                            transformation=transformation,
                            return_offset=true)
        @test display(p) != NaN
        return samp, method, fit, p, offset
    end
end
    
function partest(parname,paroffsetfact)
    samp, method, fit, p, offset = predictest()
    drift = 0.0
    down = 0.0
    for paroffset in paroffsetfact .* [-1,1]
        if parname=="drift"
            drift = fit.drift[1] + paroffset
            down = fit.down[2]
        elseif parname=="down"
            drift = fit.drift[1]
            down = fit.down[2] + paroffset
        end
        adjusted_fit = Gfit(method;
                            blank=fit.blank,
                            drift=[drift,0.0],
                            down=[0.0,down])
        plotFitted!(p,samp,method,adjusted_fit;
                    den=method.fractionation.channels.D,
                    transformation="log",
                    offset=offset,
                    linecolor="red")
    end
    @test display(p) != NaN
end

function driftest()
    partest("drift",0.5)
end

function downtest()
    partest("down",1.0)
end

function processsettings(option="Lu-Hf")
    dname = joinpath("data",option)
    head2name = true
    refmats = nothing
    snum = 1
    if option == "Lu-Hf"
        refmats =  Dict("BP" => "BP")
    elseif option == "Rb-Sr"
        refmats = Dict("MDC" => "MDC -")
        snum = 4
    elseif option == "K-Ca"
        refmats = Dict("MDC" => "MDC_")
        snum = 1 # 1, 2, 5, 6
    elseif option == "U-Pb"
        head2name = false
        refmats = Dict("Plesovice" => "STDCZ",
                       "91500" => "91500")
        snum = 3
    end
    method = getmethod(option,refmats)
    return (dname, head2name, method, refmats, snum)
end

function processtest(option="Lu-Hf";
                     show=true,
                     verbose=false,
                     transformation="log")
    dname, head2name, method, refmats, snum = processsettings(option)
    myrun = load(dname;format="Agilent",head2name=head2name)
    setGroup!(myrun,refmats)
    fit = process!(myrun,method;verbose=verbose)
    if verbose
        println("drift=",fit.drift,", down=",fit.down)
    end
    if show
        p = KJ.plot(myrun[snum],method;fit=fit,
                    transformation=transformation,
                    den=method.fractionation.channels.D)
        @test display(p) != NaN
    end
    return myrun, method, fit
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

function histest(option="Lu-Hf";show=true)
    myrun, method, fit = processtest(option)

    Pm = Float64[]; Dm = Float64[]; dm = Float64[]
    Pp = Float64[]; Dp = Float64[]; dp = Float64[]
    for samp in myrun
        if samp.group !== "sample"
            c = Cruncher(samp,method.fractionation,fit.blank)
            append!(Pm,c.pmb+c.bpt)
            append!(Dm,c.Dombi+c.bDot)
            append!(dm,c.bomb+c.bbot)
            p = predict(samp,method,fit)
            append!(Pp,p.P)
            append!(Dp,p.D)
            append!(dp,p.d)
        end
    end
    if show
        plot_residuals(Pm,Dm,dm,Pp,Dp,dp)
        df = DataFrame(Pm=Pm,Dm=Dm,dm=dm,Pp=Pp,Dp=Dp,dp=dp)
        CSV.write("output/hist.csv",df)
    end
end

function PAtest()
    dname, head2name, method, snum = processsettings("Lu-Hf")
    myrun = load(dname)
    method.PAcutoff = 1e7
    fit = process!(myrun,method)
    return myrun, method, fit
end

function atomictest(option="Lu-Hf")
    myrun, method, fit = processtest(option;show=false)
    a = atomic(myrun[2],method,fit)
    p1 = Plots.scatter(a.P,a.D,label="",xlabel="P",ylabel="D")
    p2 = Plots.scatter(a.D,a.d,label="",xlabel="D",ylabel="d")
    p = Plots.plot(p1,p2,layout=(1,2))
    @test display(p) != NaN
end

function averatest(option="Lu-Hf",verbose=false)
    myrun, method, fit = processtest(option;show=false)
    ratios = averat(myrun,method,fit)
    if verbose println(ratios) end
    return ratios, method
end

function exporttest()
    prefixes = Dict("Lu-Hf" => "hogsbo",
                    "Rb-Sr" => "EntireCreek",
                    "K-Ca" => "EntCrk",
                    "U-Pb" => "GJ1")
    for option in ["Lu-Hf","Rb-Sr","K-Ca","U-Pb"]
        ratios, method = averatest(option)
        selection = prefix2subset(ratios,prefixes[option])
        CSV.write(joinpath("output",option * ".csv"),selection)
        export2IsoplotR(selection,method,fname = joinpath("output",option * ".json"))
    end
end

function iCaptest(verbose=true)
    myrun = load("data/iCap",format="ThermoFisher")
    if verbose summarise(myrun;verbose=true,n=5) end
end

function carbonatetest(verbose=false)
    myrun = load("data/carbonate",format="Agilent")
    standards = Dict("WC1"=>"WC1")
    method = getmethod("U-Pb",standards)
    fit = process!(myrun,method)
    export2IsoplotR(myrun,method,fit;
                    prefix="Duff",
                    fname="output/Duff.json")
    p = KJ.plot(myrun[4],method;fit=fit,
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
    @test internal[2] == 1.476e6
end

function concentrationtest()
    myrun = load("data/Lu-Hf",format="Agilent")
    method = Cmethod(myrun;
                     standards=["NIST612"],
                     internal=("Al27 -> 27",1.2e5))
    setGroup!(myrun,Dict("NIST612" => "NIST612p"))
    fit = process!(myrun,method)
    conc = concentrations(myrun,method,fit)
    p = KJ.plot(myrun[4],method;fit=fit,
                transformation="log",
                den=method.internal[1])
    @test display(p) != NaN
    CSV.write("output/concentrations.csv",conc)
    return conc
end

function internochrontest(show=true)
    myrun, method, fit = processtest("Lu-Hf";show=false)
    isochron = internochron(myrun,method,fit)
    CSV.write("output/isochron.csv",isochron)
    if show
        p = internoplot(myrun[2],method,fit;xlim=[0,100])
        @test display(p) != NaN
    end
end

function internochronUPbtest(show=true)
    myrun = load("data/carbonate",format="Agilent")
    setGroup!(myrun,Dict("WC1"=>"WC1"))
    method = Gmethod("U-Pb";standards=["WC1"])
    fit = process!(myrun,method)
    isochron = internochron(myrun,method,fit)
    CSV.write("output/isochronUPb.csv",isochron)
    if show
        p = internoplot(myrun[5],method,fit)
        @test display(p) != NaN
    end
end

function maptest()
    myrun = load("data/timestamp/NHM_cropped.csv",
                 "data/timestamp/NHM_timestamps.csv";
                 format="Agilent")
    method = Cmethod(myrun;
                     standards=["NIST612"],
                     internal=getInternal("zircon","Si29"))
    setGroup!(myrun,Dict("NIST612" => "NIST612"))
    fit = process!(myrun,method)
    conc = concentrations(myrun[10],method,fit)
    p = plotMap(conc,"ppm[U] from U238";
                clims=(0,500),
                colorbar_scale=:identity)
    CSV.write("output/Umap.csv",conc)
    @test display(p) != NaN
end

function map_dating_test()
    myrun = load("data/timestamp/NHM_cropped.csv",
                 "data/timestamp/NHM_timestamps.csv";
                 format="Agilent")
    method = Gmethod("U-Pb";standards=["91500"])
    setGroup!(myrun,Dict("91500"=>"91500"))
    fit = process!(myrun,method)
    a = atomic(myrun[10],method,fit;add_xy=true)
    df = DataFrame(a)
    p = plotMap(df,"D";clims=(0,1e4))
    @test display(p) != NaN
end

function map_fail_test()
    myrun, method, fit = processtest()
    P, D, d, x, y = atomic(myrun[2],method,fit)
    df = DataFrame(P=P,D=D,d=d,x=x,y=y)
    p = plotMap(df,"P")
    @test display(p) != NaN
    conc = concentrationtest()
    p = plotMap(conc,"Lu175 -> Lu175")
    @test display(p) != NaN
end

function glass_only_test()
    myrun = load("data/U-Pb",format="Agilent",head2name=false)
    method = Gmethod("U-Pb";standards=["NIST610","NIST612"])
    setGroup!(myrun,Dict("NIST610" => "610", "NIST612" => "612"))
    fit = process!(myrun,method)
    export2IsoplotR(myrun,method,fit;
                    fname="output/UPb_with_glass.json")
    p = KJ.plot(myrun[9],method;fit=fit,
                transformation="log",den="Pb206")
    display(p)
end

function synthetictest(;drift=[0.0],down=[0.0,0.0],kw...)
    method = getmethod("Lu-Hf",Dict("BP" => "BP"))
    method.ndrift = length(drift)
    method.ndown = length(down)-1
    myrun, fit = synthetic!(method;
                            drift=drift,
                            down=down,
                            lambda=1.867e-5,
                            t_std=1745.0,
                            t_smp=1029.7,
                            y0_smp=3.55,
                            y0_glass=3.544842,
                            kw...)
    return myrun, method, fit
end

function SS4test(run::Vector{Sample},
                 method::Gmethod,
                 fit::Gfit)
    out = 0.0
    for samp in run
        if samp.group in method.fractionation.standards
            a = getAnchor(method.name,samp.group)
            c = Cruncher(samp,method.fractionation,fit.blank)
            ft = polyFac(fit.drift,c.t)
            FT = polyFac(fit.down,c.T)
            out += SS(a,ft,FT;c...)
        end
    end
    return out
end

function SStest()
    myrun, method, truefit = synthetictest(;drift=[0.0],down=[0.0,0.0])
    nstep = 50
    driftss = fill(0.0,nstep)
    driftfit = deepcopy(truefit)
    dd = range(start=truefit.drift[1]-1.0,stop=truefit.drift[1]+1.0,length=nstep)
    for i in eachindex(dd)
        driftfit.drift[1] = driftfit.adrift[1] = dd[i]
        driftss[i] = SS4test(myrun,method,driftfit)
    end
    p = Plots.plot(dd,driftss,seriestype=:line,label="drift")
    dwn = range(start=truefit.down[2]-1.0,stop=truefit.down[2]+1.0,length=nstep)
    downss = fill(0.0,nstep)
    downfit = deepcopy(truefit)
    for i in eachindex(dwn)
        downfit.down[2] = dwn[i]
        downss[i] = SS4test(myrun,method,downfit)
    end
    Plots.plot!(p,dwn,downss,seriestype=:line,linecolor=:red,label="down")
    display(p)
end

function accuracytest(;drift=[0.0],down=[0.0,0.0],show=true,kw...)
    myrun, method, truefit = synthetictest(;drift=drift,down=down,kw...)
    fit = process!(myrun,method)
    lldrift = fit.drift[1] - 3*sqrt(fit.covmat[1,1])
    uldrift = fit.drift[1] + 3*sqrt(fit.covmat[1,1])
    @test uldrift > truefit.drift[1] && lldrift < truefit.drift[1]
    lldown = fit.down[2] - 3*sqrt(fit.covmat[2,2])
    uldown = fit.down[2] + 3*sqrt(fit.covmat[2,2])
    @test uldown > truefit.down[2] && lldown < truefit.down[2]
    if show
        den = nothing # "Hf176 -> 258" #
        p1 = KJ.plot(myrun[1],method;fit=fit,
                     transformation="sqrt",den=den)
        p2 = KJ.plot(myrun[4],method;fit=fit,
                     transformation="sqrt",den=den)
        p = Plots.plot(p1,p2,layout=(1,2))
        @test display(p) != NaN
    end
end

function biastest()
    method = Gmethod("Re-Os";
                     ions = (P="Re187",D="Os187",d="Os188"),
                     proxies = (P="Re185",D="Os187",d="Os189"),
                     channels = (P="Re185 -> 185",D="Os187 -> 251",d="Os189 -> 253"),
                     standards = ["QMolyHill"],
                     bias = Dict("Os" => ["NiS-3"]),
                     nbias = 2)
    method.interference = Interference(;ions = Dict("Os187" => ["Re187"]),
                                        proxies = Dict("Re187" => "Re185"),
                                        channels = Dict("Re185" => "Re185 -> 249"),
                                        bias = Dict("Re" => ["NIST612"]))
    myrun = load("data/Re-Os";format="Agilent")
    refmats = Dict("QMolyHill" => "Qmoly",
                   "NiS-3" => "Nis3",
                   "NIST612" => "Nist_massbias")
    setGroup!(myrun,refmats)
    blanks = fitBlanks(myrun)
    bias_interference = interference_bias(myrun,method,blanks)
    bias_fractionation = fractionation_bias(myrun,method,blanks)
end

function interference_test()
    
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

function dependencytest()
    Aqua.test_all(KJ;undocumented_names=false)
end

Plots.closeall()

# @testset "load" begin loadtest(;verbose=true) end
# @testset "plot raw data" begin plottest(2) end
# @testset "set selection window" begin windowtest() end
# @testset "set method and blanks" begin blanktest() end
# @testset "moving median test" begin mmediantest() end
# @testset "outlier detection" begin outliertest_synthetic() end
# @testset "outlier detection" begin outliertest_sample() end
# @testset "create method" begin methodtest() end
# @testset "assign refmats" begin refmattest(true) end
# @testset "predict Lu-Hf" begin predictest("Lu-Hf";snum=1) end
# @testset "predict Rb-Sr" begin predictest("Rb-Sr";snum=2) end
# @testset "predict K-Ca" begin predictest("K-Ca";snum=1) end
# @testset "predict drift" begin driftest() end
# @testset "predict down" begin downtest() end
# @testset "Lu-Hf" begin processtest("Lu-Hf") end
# @testset "Rb-Sr" begin processtest("Rb-Sr") end
# @testset "K-Ca" begin processtest("K-Ca") end
# @testset "U-Pb" begin processtest("U-Pb") end
# @testset "hist" begin histest() end
# @testset "PA test" begin PAtest() end
# @testset "atomic test" begin atomictest("Rb-Sr") end
# @testset "averat test" begin averatest("K-Ca") end
# @testset "export" begin exporttest() end
# @testset "iCap" begin iCaptest() end
# @testset "carbonate" begin carbonatetest() end
# @testset "timestamp" begin timestamptest() end
# @testset "stoichiometry" begin mineraltest() end
# @testset "concentration" begin concentrationtest() end
# @testset "Lu-Hf internochron" begin internochrontest() end
# @testset "UPb internochron" begin internochronUPbtest() end
# @testset "concentration map" begin maptest() end
# @testset "isotope ratio map" begin map_dating_test() end
# @testset "map fail test" begin map_fail_test() end
# @testset "glass as age standard test" begin glass_only_test() end
# @testset "extension test" begin extensiontest() end
# @testset "synthetic data" begin SStest() end
# @testset "accuracy test 1" begin accuracytest() end
# @testset "accuracy test 2" begin accuracytest(drift=[-2.0]) end
# @testset "accuracy test 3" begin accuracytest(down=[0.0,0.5]) end
@testset "bias test" begin biastest() end
# @testset "interference test" begin interference_test() end
# @testset "TUI test" begin TUItest() end
# @testset "dependency test" begin dependencytest() end

# TUI()