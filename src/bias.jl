function fit_bias(run::Vector{Sample},
                  method::Gmethod,
                  element::AbstractString,
                  blank::AbstractDataFrame)

    calibration = method.bias[element]
    cruncher_groups = Dict()
    for group in calibration.standards
        standard = method.groups[group]
        y = getAnchor(method.name,standard).y
        selection = group2selection(run,group)
        ns = length(selection)
        crunchers = Vector{NamedTuple}(undef,ns)
        for i in eachindex(selection)
            crunchers[i] = BCruncher(run[i],method,element,blank)
        end
        cruncher_groups[standard] = (y=y,crunchers=crunchers)
    end

    m1 = get_proxy_isotope(calibration.num.ion,element)
    m2 = get_proxy_isotope(calibration.den.ion,element)

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
                   element::AbstractString,
                   blank::AbstractDataFrame)

    dat = swinData(samp)
    calibration = method.bias[element]
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

    Ib = fill(1.0,length(t))
    ID = fill(1.0,length(t))

    bd = iratio(method.d.proxy,method.d.ion)

    return (bmb=bmb,Dmb=Dmb,vb=vb,vD=vD,sDb=sDb,bd=bd,ID=ID,Ib=Ib,t=t)

end