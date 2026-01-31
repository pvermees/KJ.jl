function fit_bias(run::Vector{Sample},
                  method::Gmethod,
                  blank::AbstractDataFrame)
    out = DataFrame()
    pairings = (method.P,method.D,method.d)
    for pairing in pairings
        for interference in pairing.interferences
            calibration = interference.bias
            if length(calibration.standards) > 0
                bias_key = channel2element(calibration.num.ion)
                out[!,bias_key] = fit_bias(run,method,calibration,blank)
            end
        end
    end
    calibration = method.bias
    if length(calibration.standards) > 0
        bias_key = channel2element(calibration.num.ion)
        out[!,bias_key] = fit_bias(run,method,calibration,blank;
                                   interference_bias=out)
    end
    return out
end
function fit_bias(run::Vector{Sample},
                  method::Gmethod,
                  calibration::Calibration,
                  blank::AbstractDataFrame;
                  interference_bias::AbstractDataFrame=DataFrame())

    cruncher_groups = Dict()
    for group in calibration.standards
        standard = method.groups[group]
        y = iratio(calibration.num.ion,calibration.den.ion)
        if isnothing(y)
            y = getAnchor(method.name,standard).y
        end
        selection = group2selection(run,group)
        ns = length(selection)
        crunchers = Vector{NamedTuple}(undef,ns)
        for i in eachindex(selection)
            crunchers[i] = BCruncher(run[i],method,calibration,blank;
                                     interference_bias=interference_bias)
        end
        cruncher_groups[standard] = (y=y,crunchers=crunchers)
    end

    m1 = get_proxy_isotope(calibration.num.ion)
    m2 = get_proxy_isotope(calibration.den.ion)

    init = fill(0.0,method.nbias)
    objective = (par) -> SS(par,m1,m2,cruncher_groups)
    optimum = Optim.optimize(objective,init)

    return Optim.minimizer(optimum)
    
end
export fit_bias

function bias_correction(par::AbstractVector,
                         mass1::Int,
                         mass2::Int;
                         t::AbstractVector,
                         other...)
    beta = (1/mass1 - 1/(mass1+1))/(1/mass1 - 1/mass2)
    return polyFac(par,t).^beta
end

function BCruncher(samp::Sample,
                   method::Gmethod,
                   calibration::Calibration,
                   blank::AbstractDataFrame;
                   interference_bias::AbstractDataFrame=DataFrame())

    dat = swinData(samp)
    bm = dat[:,calibration.num.channel]
    Dm = dat[:,calibration.den.channel]

    t = dat.T
    bbt = polyVal(blank[:,method.d.channel],t)
    bDt = polyVal(blank[:,method.D.channel],t)

    Dmb = Dm - bDt
    bmb = bm - bbt

    sig = hcat(Dmb,bmb)
    covmat = df2cov(sig)
    vD = covmat[1,1]
    vb = covmat[2,2]
    sDb = covmat[1,2]

    bd = iratio(method.d.proxy,method.d.ion)
    Ib = interference_correction(dat,method.d.interferences;
                                 bias=interference_bias,blank=blank)
    ID = interference_correction(dat,method.D.interferences;
                                 bias=interference_bias,blank=blank)

    return (bmb=bmb,Dmb=Dmb,vb=vb,vD=vD,sDb=sDb,bd=bd,ID=ID,Ib=Ib,t=t)

end