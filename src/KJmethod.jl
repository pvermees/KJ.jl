function default_ions(name)
    m = get(_KJ["methods"],name)
    return (P=String(m.P),D=String(m.D),d=String(m.d))
end

function channel2proxy(channel::AbstractString)
    all_elements = string.(keys(_KJ["nuclides"]))
    matching_elements = filter(x -> occursin(x, channel), all_elements)
    matching_element = argmax(length,matching_elements)
    matching_isotope = get_proxy_isotope(channel;element=matching_element)
    if length(matching_isotope) > 0
        return matching_element * string(matching_isotope)
    else
        return nothing
    end
end

function get_proxy_isotope(channel::AbstractString;
                           element::AbstractString=channel2element(channel))
    all_isotopes = _KJ["nuclides"][element]
    matching_isotope = filter(x -> occursin(string(x), channel), all_isotopes)
    return matching_isotope[1]
end

function Cmethod(run::Vector{Sample};
                 groups::AbstractDict=Dict{String,String}(),
                 internal::Tuple=(nothing,nothing),
                 nblank::Int=2)
    ch = getChannels(run)
    el = channel2element.(ch)
    elements = NamedTuple{Tuple(Symbol.(ch))}(Tuple(el))
    return Cmethod(elements,groups,internal,nblank)
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