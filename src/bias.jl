function Cruncher(samp::Sample,
                  num_channel::AbstractString,
                  den_channel::AbstractString,
                  blank::AbstractDataFrame)

    dat = swinData(samp)
    
    dm = dat[:,num_channel]
    Dm = dat[:,den_channel]

    t = dat.t

    bdt = polyVal(blank[:,num_channel],t)
    bDt = polyVal(blank[:,den_channel],t)

    Dmb = Dm - bDt
    dmb = dm - bdt

    sig = hcat(Dmb,dmb)
    covmat = df2cov(sig)
    vD = covmat[1,1]
    vd = covmat[2,2]
    sDd = covmat[1,2]

    return (Dm=Dm,dm=dm,
            bDt=bDt,bdt=bdt,
            Dmb=Dmb,dmb=dmb,
            vD=vD,vd=vd,sDd=sDd,
            t=t)

end

function bias_prep(F::Fractionation)
    element = channel2element(F.ions.D)
    standards = F.bias[element]
    num = (ion=F.proxies.d,channel=F.channels.d)
    den = (ion=F.proxies.D,channel=F.channels.D)
    return (standards=standards, num=num, den=den)
end

function bias_prep(I::Interference,
                   target_channel::AbstractString,
                   interfering_ion::AbstractString)
    interfering_element = channel2element(interfering_ion)
    interference_proxy = I.proxies[interfering_ion]
    interference_proxy_channel = I.channels[interference_proxy]
    standards = I.bias[interfering_element]
    num = (ion=interference_proxy,channel=interference_proxy_channel)
    den = (ion=interfering_ion,channel=target_channel)
    return (standards=standards, num=num, den=den)
end

function fractionation_cruncher_groups(run::Vector{Sample},
                                       method::Gmethod,
                                       blank::AbstractDataFrame)
    out = Dict()
    standards, num, den = bias_prep(method.fractionation)
    element = channel2element(method.fractionation.ions.D)
    out[element] = bias_cruncher_groups_helper(run,method.name,blank,standards,num,den)
    return out
end

function bias_cruncher_groups_helper(run::Vector{Sample},
                                     methodname::AbstractString,
                                     blank::AbstractDataFrame,
                                     standards::AbstractVector,
                                     num::NamedTuple{(:ion,:channel)},
                                     den::NamedTuple{(:ion,:channel)})
    cruncher_groups = Dict()
    for standard in standards
        element = channel2element(num.ion)
        y = iratio(element,num.ion,den.ion)
        if isnothing(y)
            y = getAnchor(methodname,standard)
        end
        selection = group2selection(run,standard)
        ns = length(selection)
        crunchers = Vector{NamedTuple}(undef,ns)
        for i in eachindex(selection)
            crunchers[i] = Cruncher(run[selection[i]],
                                    num.channel,den.channel,
                                    blank)
        end
        cruncher_groups[standard] = (anchor=y,crunchers=crunchers)
    end
    return cruncher_groups
end

function interference_cruncher_groups(run::Vector{Sample},
                                      method::Gmethod,
                                      blank::AbstractDataFrame)
    out = Dict()
    F = method.fractionation
    Fchannels = Dict(zip(values(F.proxies),values(F.channels)))
    I = method.interference
    for (target,interfering_ions) in I.ions
        target_channel = Fchannels[target]
        for interfering_ion in interfering_ions
            standards, num, den = bias_prep(I,target_channel,interfering_ion)
            cruncher_groups = bias_cruncher_groups_helper(run,method.name,blank,standards,num,den)
            interfering_element = channel2element(interfering_ion)
            out[interfering_element] = cruncher_groups
        end
    end
    return out
end

function fractionation_bias(run::Vector{Sample},
                            method::Gmethod,
                            blank::AbstractDataFrame)
    cruncher_groups = fractionation_cruncher_groups(run,method,blank)
    return bias(cruncher_groups,method.nbias)
end
export fractionation_bias

function interference_bias(run::Vector{Sample},
                           method::Gmethod,
                           blank::AbstractDataFrame)
    cruncher_groups = interference_cruncher_groups(run,method,blank)
    return bias(cruncher_groups,method.nbias)
end
export interference_bias

function bias(cruncher_groups::AbstractDict,
              nbias::Int)
    out = DataFrame()
    for (element,cruncher_groups) in cruncher_groups
        init = fill(0.0,nbias)
        objective = (par) -> SS(par,cruncher_groups)
        optimum = Optim.optimize(objective,init)
        solution = Optim.minimizer(optimum)
        out[:,element] = solution
    end
    return out
end
export bias