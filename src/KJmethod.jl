function default_ions(name)
    m = get(_KJ["methods"], name)
    return (P=String(m.P), D=String(m.D), d=String(m.d))
end

function channel2proxy(channel::AbstractString;
    elements::AbstractVector=string.(keys(_KJ["nuclides"])))
    matching_elements = filter(x -> occursin(x, channel), elements)
    if length(matching_elements) < 1
        return nothing
    end
    matching_element = argmax(length, matching_elements)
    matching_isotope = get_proxy_isotope(channel; element=matching_element)
    if isnothing(matching_isotope)
        return nothing
    else
        return matching_element * string(matching_isotope)
    end
end

function get_proxy_isotope(channel::AbstractString;
    element::AbstractString=channel2element(channel))
    all_isotopes = _KJ["nuclides"][element]
    matching_isotope = filter(x -> occursin(string(x), channel), all_isotopes)
    if length(matching_isotope) > 0
        return matching_isotope[1]
    else
        return nothing
    end
end

"""
    Cmethod(run::Vector{Sample}; groups=Dict{String,String}(), internal=(nothing, nothing), nblank=2)

Build a concentration method from sample channels by inferring the element
for each channel, then delegating to the typed `Cmethod` constructor.
"""
function Cmethod(run::Vector{Sample};
    groups::AbstractDict=Dict{String,String}(),
    internal::Dict{String,Tuple{String,ℝ}}=Dict{String,Tuple{String,ℝ}}(),
    nblank::Int=2) where ℝ <: Real
    ch = getChannels(run)
    el = channel2element.(ch)
    elements = NamedTuple{Tuple(Symbol.(ch))}(Tuple(el))
    return Cmethod(elements, groups, internal, nblank)
end

function getConcentrations(method::Cmethod,
    refmat::AbstractString)
    all_concs = get(_KJ["glass"], refmat)
    els_rm = names(all_concs)
    channels = getChannels(method)
    out = DataFrame(zeros(1, length(channels)), channels)
    for (ch, el) in pairs(method.elements)
        if in(el, els_rm)
            out[1, ch] = all_concs[el]
        else
            @warn "$el missing from reference material ($refmat) concentrations. Output for this element will be NaN."
            out[1, ch] = NaN
        end
    end
    return out
end
