function Gmethod(name::String;
                 groups::AbstractDict=Dict{String,String}(),
                 P::Pairing=Pairing(ion=default_ions(name).P),
                 D::Pairing=Pairing(ion=default_ions(name).D),
                 d::Pairing=Pairing(ion=default_ions(name).d),
                 bias::AbstractDict=Dict(),
                 standards::AbstractSet=Set{String}(),
                 nblank::Int=2,
                 ndrift::Int=2,
                 ndown::Int=1,
                 nbias::Int=1,
                 PAcutoff::Float64=Inf)
    return Gmethod(name,groups,P,D,d,bias,standards,
                   nblank,ndrift,ndown,nbias,PAcutoff)
end

function default_ions(name)
    m = get(_KJ["methods"],name)
    return (P=String(m.P),D=String(m.D),d=String(m.d))
end

function channels2proxies(channels::NamedTuple{(:P,:D,:d)})
    return (P=channel2proxy(channels.P),
            D=channel2proxy(channels.D),
            d=channel2proxy(channels.d))
end
function channel2proxy(channel::AbstractString)
    proxy = nothing
    all_elements = string.(keys(_KJ["nuclides"]))
    matching_elements = filter(x -> occursin(x, channel), all_elements)
    matching_element = argmax(length,matching_elements)
    proxy = get_proxy_isotope(channel,matching_element)
    return proxy
end

function get_proxy_isotope(channel::AbstractString,
                           matching_element::AbstractString)
    all_isotopes = string.(_KJ["nuclides"][matching_element])
    matching_isotope = filter(x -> occursin(x, channel), all_isotopes)
    if length(matching_isotope) == 1
        return matching_element * matching_isotope[1]
    else
        return nothing
    end
end

function Cmethod(run::Vector{Sample};
                 standards::AbstractVector=String[],
                 internal::Tuple=(nothing,nothing),
                 nblank::Int=2)
    ch = getChannels(run)
    el = channel2element.(ch)
    elements = NamedTuple{Tuple(Symbol.(ch))}(Tuple(el))
    return Cmethod(elements,Set(standards),internal,nblank)
end

function getConcentrations(method::Cmethod,
                           refmat::AbstractString)
    all_concs = get(_KJ["glass"],refmat)
    channels = getChannels(method)
    out = DataFrame(zeros(1, length(channels)), channels)
    for (ch,el) in pairs(method.elements)
        out[1,ch] = all_concs[el]
    end
    return out
end