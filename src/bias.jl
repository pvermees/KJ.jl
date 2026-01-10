function bias(run::Vector{Sample},
              method::Gmethod,
              blank::AbstractDataFrame)
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
            cruncher_groups = Dict()
            for standard in standards
                y = iratio(interfering_element,
                           interfering_ion,
                           interference_proxy)
                if isnothing(y)
                    y = getAnchor(method.name,standard)
                end
                selection = group2selection(run,standard)
                ns = length(selection)
                crunchers = Vector{NamedTuple}(undef,ns)
                 for i in eachindex(selection)
                    crunchers[i] = Cruncher(run[selection[i]],
                                            target_channel,
                                            interference_proxy_channel,
                                            blank)
                end
                cruncher_groups[standard] = (anchor=y,crunchers=crunchers)
            end
        end
    end
    return nothing
end

function bias(run::Vector{Sample},
              fractionation::Fractionation)
    c = Cruncher(run,fractionation)
     return nothing
end
export bias

function Cruncher(samp::Sample,
                  interference_channel::AbstractString,
                  proxy_channel::AbstractString,
                  blank::AbstractDataFrame)

    dat = swinData(samp)
    
    Dm = dat[:,interference_channel]
    dm = dat[:,proxy_channel]

    t = dat.t

    bDt = polyVal(blank[:,interference_channel],t)
    bdt = polyVal(blank[:,proxy_channel],t)

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
            vD=vD,vd=vd,sDd=sDd)

end