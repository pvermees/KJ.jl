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
                selection = group2selection(run,standard)
                ns = length(selection)
                crunchers = Vector{BCruncher}(undef,ns)
            end
        end
    end
    return nothing
end

function bias(run::Vector{Sample},
              fractionation::Fractionation)
    c = BCruncher(run,fractionation)
     return nothing
end
export bias

function BCruncher(run::Vector{Sample},
                   method::Gmethod,
                   blank::AbstractDataFrame)
    return nothing
end

function BCruncher(run::Vector{Sample},
                   fractionation::Fractionation)
    return nothing
end