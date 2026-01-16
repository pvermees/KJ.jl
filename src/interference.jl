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
                                  method::KJmethod)
    for (target,interferences) in method.interferences
        i = findfirst(==(target),unpack(method.proxies))
        target_channel = unpack(method.channels)[i]
        for interference in interferences
            ratio = iratio(interference.ion,interference.proxy)
            correction = ratio .* samp.dat[:,interference.channel]
            samp.dat[:,target_channel] .-= correction
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