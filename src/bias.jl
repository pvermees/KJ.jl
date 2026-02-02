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
        out[bias_key] = fit_bias(run,method,calibration,blank)
    end
    return out
end
function fit_bias(run::Vector{Sample},
                  m::Gmethod,
                  c::Calibration,
                  blank::AbstractDataFrame)

    mass_num = get_proxy_isotope(c.num.ion)
    mass_den = get_proxy_isotope(c.den.ion)
    bd = calibration2bd(m,c)

    cruncher_groups = Dict()
    for group in c.standards
        standard = m.groups[group]
        y = get_bias_truth(m,c,standard)
        selection = group2selection(run,group)
        ns = length(selection)
        crunchers = Vector{NamedTuple}(undef,ns)
        for i in eachindex(selection)
            samp = run[selection[i]]
            crunchers[i] = BCruncher(samp,c,blank)
        end
        cruncher_groups[standard] = (y=y,crunchers=crunchers)
    end

    init = fill(0.0,m.nbias)
    objective = (par) -> SS(par,mass_num,mass_den,bd,cruncher_groups)
    optimum = Optim.optimize(objective,init)
    fit = Optim.minimizer(optimum)

    return Bias(mass_num,mass_den,fit)
end
export fit_bias

function calibration2bd(m::Gmethod,
                        c::Calibration)
    if c.den.ion == m.D.ion && c.num.ion != m.d.ion
        return iratio(c.num.ion,m.d.ion)
    else
        return 1.0
    end
end

function get_bias_truth(m::Gmethod,
                        c::Calibration,
                        standard::AbstractString)
    y = iratio(c.num.ion,c.den.ion)
    if isnothing(y)
        y = getAnchor(m.name,standard).y
    end
    return y
end

function bias_correction(bias::Bias,
                         mass_num::Int,
                         mass_den::Int;
                         t::AbstractVector,
                         other...)
    beta = log(mass_num/mass_den)/log(bias.mass_num/bias.mass_den)
    return polyFac(bias.par,t).^beta
end
function bias_correction(bias::Bias,
                         num::AbstractString,
                         den::AbstractString,
                         t::AbstractVector)
    mass_num = get_proxy_isotope(num)
    mass_den = get_proxy_isotope(den)
    return bias_correction(bias,mass_num,mass_den;t=t)
end

function BCruncher(samp::Sample,
                   calibration::Calibration,
                   blank::AbstractDataFrame)

    dat = swinData(samp)
    t = dat.T
    Dch = calibration.den.channel
    bch = calibration.num.channel

    Dm = dat[:,Dch]
    bm = dat[:,bch]
    bDt = polyVal(blank[:,Dch],t)
    bbt = polyVal(blank[:,bch],t)
    Dmb = Dm - bDt
    bmb = bm - bbt
    b_vs_D = hcat(Dmb,bmb)
    covmat = df2cov(b_vs_D)
    vD = covmat[1,1]
    vb = covmat[2,2]
    sDb = covmat[1,2]

    return (Dmb=Dmb,bmb=bmb,
            bDt=bDt,bbt=bbt,
            vD=vD,vb=vb,sDb=sDb,t=t)

end