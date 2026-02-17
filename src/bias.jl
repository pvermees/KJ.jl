function fit_bias(run::Vector{Sample},
                  method::Gmethod,
                  blank::AbstractDataFrame)
    out = Dict{String,AbstractBias}()
    pairings = (method.P,method.D,method.d)
    for pairing in pairings
        for (key,interference) in pairing.interferences
            add_bias!(out,run,method,blank,key,interference)
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

    init = [init_bias(cruncher_groups);fill(0.0,m.ndrift-1)]
    objective = (par) -> SS(par,mass_num,mass_den,bd,cruncher_groups)
    optimum = Optim.optimize(objective,init)
    fit = Optim.minimizer(optimum)
    return Bias(mass_num,mass_den,fit)
end
export fit_bias

function init_bias(cruncher_groups::AbstractDict)
    bias = 0.0
    ncg = length(cruncher_groups)
    for group in values(cruncher_groups)
        ytrue = group.y
        bmb = 0.0
        Dmb = 0.0
        for c in group.crunchers
            bmb += sum(c.bmb)
            Dmb += sum(c.Dmb)
        end
        ymeas = bmb/Dmb
        bias += log(ymeas/ytrue)/ncg
    end
    return bias
end

function add_bias!(bias::AbstractDict,
                   run::Vector{Sample},
                   method::Gmethod,
                   blank::AbstractDataFrame,
                   ion::AbstractString,
                   interference::Interference)
    calibration = interference.bias
    if length(calibration.standards) > 0
        bias_key = channel2element(calibration.num.ion)
        bias[bias_key] = fit_bias(run,method,calibration,blank)
    end
end

function calibration2bd(m::Gmethod,
                        c::Calibration)
    if c.den.ion == m.D.ion && c.num.ion != m.d.ion
        return iratio(c.num.ion,m.d.ion)
    elseif c.num.ion == m.D.ion && c.den.ion != m.d.ion
        return iratio(m.d.ion,c.den.ion) # flipped
    else
        return 1.0
    end
end

function get_bias_truth(m::Gmethod,
                        c::Calibration,
                        standard::AbstractString)
    y = iratio(c.num.ion,c.den.ion)
    if isnothing(y)
        a = getAnchor(m.name,standard)
        if m.D.proxy == c.den.ion
            y = a.y
        elseif m.D.proxy == c.num.ion
            y = 1/a.y
        else
            @warn "The calibration isotope is not part of the method."
            y = 1.0
        end
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
    Dch = calibration.den.channel
    bch = calibration.num.channel
    return BCruncher(samp,Dch,bch,blank)
end

function BCruncher(samp::Sample,
                   Dch::AbstractString,
                   bch::AbstractString,
                   blank::AbstractDataFrame)
    dat = swinData(samp)
    t = dat.T

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