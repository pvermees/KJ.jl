function Interference(;ions::Dict{String,Vector{String}}=Dict{String,Vector{String}}(),
                       channels::Dict{String,String}=Dict{String,String}(),
                       proxies::Dict{String,String}=Dict{String,String}(),
                       bias::Dict{String,Vector{String}}=Dict{String,Vector{String}}())
    return Interference(ions,channels,proxies,bias)
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