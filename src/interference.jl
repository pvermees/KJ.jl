function Interference(;ion::AbstractString="",
                       proxy::AbstractString="",
                       channel::AbstractString="",
                       bias_key::AbstractString=channel2element(ion))
    return Interference(ion,proxy,channel,bias_key)
end

function REEInterference(;proxychannel::AbstractString="",
                          numchannel::AbstractString="",
                          denchannel::AbstractString="")
    return REEInterference(proxychannel,numchannel,denchannel)
end
    
function interference_correction!(run::Vector{Sample},
                                  method::KJmethod;
                                  bias::AbstractDataFrame=init_bias(method))
    for samp in run
        interference_correction!(samp,method;bias=bias)
    end
end
export interference_correction!

function interference_correction!(samp::Sample,
                                  method::Gmethod;
                                  bias::AbstractDataFrame=init_bias(method))
    F = method.fractionation
    I = method.interference
    for (key,proxy) in pairs(F.proxies)
        if proxy in keys(I.ions)
            target_channel = F.channels[key]
            interferences = I.ions[proxy]
            for interference_ion in interferences
                interference_element = channel2element(interference_ion)
                interference_proxy = I.proxies[interference_ion]
                interference_channel = I.channels[interference_proxy]
                mf = polyFac(bias[:,interference_element],samp.dat.t)
                ratio = mf .* iratio(interference_ion,interference_proxy)
                correction = ratio .* samp.dat[:,interference_channel]
                samp.dat[:,target_channel] .-= correction
            end
        end
    end
end

function infer_interference(method::Gmethod,
                            standard::AbstractString)
    F = method.fractionation
    I = method.interference
    Fchannels = Dict(zip(values(F.proxies),values(F.channels)))
    for (target,interfering_ions) in I.ions
        target_channel = Fchannels[target]
        for interfering_ion in interfering_ions
            proxy = I.proxies[interfering_ion]
            element = channel2element(proxy)
            if standard in I.bias[element]
                return target_channel, interfering_ion
            end
        end
    end
    return nothing, nothing
end