function bias(run::Vector{Sample},
              method::Gmethod,
              blank::AbstractDataFrame)
    out = Dict()
    F = method.fractionation
    Fchannels = Dict(zip(values(F.proxies),values(F.channels)))
    I = method.interference
    for (target,interfering_ions) in I.ions
        target_channel = Fchannels[target]
        for interfering_ion in interfering_ions
            interfering_element = channel2element(interfering_ion)
            interference_proxy = I.proxies[interfering_ion]
            interference_proxy_channel = I.channels[interference_proxy]
            standards = I.bias[interfering_element]
            num = (ion=interference_proxy,channel=interference_proxy_channel)
            den = (ion=interfering_ion,channel=target_channel)
            cruncher_groups = get_bias_cruncher_groups(run,method.name,blank,standards,num,den)
            out[interfering_element] = bias(cruncher_groups,method.nbias)
        end
    end
    return out
end

function bias(run::Vector{Sample},
              fractionation::Fractionation,
              blank::AbstractDataFrame)
    out = Dict()
    element = channel2element(fractionation.ions.D)
    standards = fractionation.bias[element]
    for standard in standards
        @infiltrate
    end
    return out
end

function bias(cruncher_groups::AbstractDict,
              nbias::Int)
    init = fill(0.0,nbias)
    objective = (par) -> SS(par,cruncher_groups)
    optimum = Optim.optimize(objective,init)
    solution = Optim.minimizer(optimum)
    return solution
end
export bias

function get_bias_cruncher_groups(run::Vector{Sample},
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