function Interference(;ions::Dict{String,Vector{String}}=Dict{String,Vector{String}}(),
                       proxies::Dict{String,String}=Dict{String,String}(),
                       channels::Dict{String,String}=Dict{String,String}(),
                       bias::Dict{String,Vector{String}}=Dict{String,Vector{String}}())
    return Interference(ions,proxies,channels,bias)
end
    
function interference_correction!(run::Vector{Sample},
                                  method::KJmethod)
    for samp in run
        interference_correction!(samp,method)
    end
end
export interference_correction!

function interference_correction!(samp::Sample,
                                  method::Gmethod)
    F = method.fractionation
    I = method.interference
    for (key,proxy) in pairs(F.proxies)
        if proxy in keys(I.ions)
            target_channel = F.channels[key]
            interferences = I.ions[proxy]
            for interference_ion in interferences
                interference_proxy = I.proxies[interference_ion]
                interference_channel = I.channels[interference_proxy]
                bias = 1.0 # TODO
                ratio = bias .* iratio(interference_ion,interference_proxy)
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