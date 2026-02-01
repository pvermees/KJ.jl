function fit_bias(run::Vector{Sample},
                  method::Gmethod,
                  blank::AbstractDataFrame)
    out = Dict{String,Bias}()
    pairings = (method.P,method.D,method.d)
    for pairing in pairings
        for interference in pairing.interferences
            calibration = interference.bias
            if length(calibration.standards) > 0
                bias_key = channel2element(calibration.num.ion)
                out[bias_key] = fit_bias(run,method,calibration,blank)
            end
        end
    end
    calibration = method.bias
    if length(calibration.standards) > 0
        bias_key = channel2element(calibration.num.ion)
        out[bias_key] = fit_bias(run,method,calibration,blank;
                                 interference_bias=out)
    end
    return out
end
function fit_bias(run::Vector{Sample},
                  method::Gmethod,
                  calibration::Calibration,
                  blank::AbstractDataFrame;
                  interference_bias::AbstractDict=Dict{String,Bias}())

    cruncher_groups = Dict()
    for group in calibration.standards
        standard = method.groups[group]
        y, m1, m2 = get_bias_truth(method,calibration,standard)
        selection = group2selection(run,group)
        ns = length(selection)
        crunchers = Vector{NamedTuple}(undef,ns)
        for i in eachindex(selection)
            crunchers[i] = BCruncher(run[i],method,calibration,blank;
                                     interference_bias=interference_bias)
        end
        cruncher_groups[standard] = (y=y,m1=m1,m2=m2,crunchers=crunchers)
    end

    m1 = get_proxy_isotope(calibration.den.ion)
    m2 = get_proxy_isotope(calibration.num.ion)

    init = fill(0.0,method.nbias)
    objective = (par) -> SS(par,m1,m2,cruncher_groups)
    optimum = Optim.optimize(objective,init)

    fit = Optim.minimizer(optimum)
    return Bias(m1,m2,fit)
end
export fit_bias

function get_bias_truth(method::Gmethod,
                        calibration::Calibration,
                        standard::AbstractString)
    y = iratio(calibration.num.ion,calibration.den.ion)
    if isnothing(y)
        y = getAnchor(method.name,standard).y
        m1 = get_proxy_isotope(method.D.ion)
        m2 = get_proxy_isotope(method.d.ion)
    else
        m1 = get_proxy_isotope(calibration.den.ion)
        m2 = get_proxy_isotope(calibration.num.ion)
    end
    return y, m1, m2
end

function bias_correction(bias::Bias,
                         m1::Int,
                         m2::Int;
                         t::AbstractVector,
                         other...)
    beta = log(m2/m1)/log(bias.m2/bias.m1)
    return polyFac(bias.par,t).^beta
end
function bias_correction(bias::Bias,
                         num::AbstractString,
                         den::AbstractString,
                         t::AbstractVector)
    m1 = get_proxy_isotope(num)
    m2 = get_proxy_isotope(den)
    return bias_correction(bias,m1,m2;t=t)
end

function BCruncher(samp::Sample,
                   method::Gmethod,
                   calibration::Calibration,
                   blank::AbstractDataFrame;
                   interference_bias::AbstractDict=Dict{String,Bias}())

    dat = swinData(samp)
    bm = dat[:,calibration.num.channel]
    Dm = dat[:,calibration.den.channel]

    t = dat.T
    bbt = polyVal(blank[:,method.d.channel],t)
    bDt = polyVal(blank[:,method.D.channel],t)

    Dmb = Dm - bDt
    bmb = bm - bbt

    b_vs_D = hcat(Dmb,bmb)
    covmat = df2cov(b_vs_D)
    vD = covmat[1,1]
    vb = covmat[2,2]
    sDb = covmat[1,2]

    bd = iratio(method.d.proxy,method.d.ion)
    Ib = interference_correction(dat,method.d.interferences;
                                 bias=interference_bias,blank=blank)
    ID = interference_correction(dat,method.D.interferences;
                                 bias=interference_bias,blank=blank)

    return (bmb=bmb,Dmb=Dmb,
            bbt=bbt,bDt=bDt,
            vb=vb,vD=vD,sDb=sDb,
            bd=bd,ID=ID,Ib=Ib,t=t)

end