function Gmethod(name::String,
                 standards::AbstractDict;
                 ions::NamedTuple{(:P,:D,:d)}=default_ions(name),
                 proxies::NamedTuple{(:P,:D,:d)}=ions,
                 channels::NamedTuple{(:P,:D,:d)}=proxies,
                 interferences::AbstractDict=Dict(),
                 nblank::Int=2,
                 ndrift::Int=2,
                 ndown::Int=1,
                 PAcutoff::Union{Nothing,Float64}=nothing,
                 glass::AbstractDict=Dict())
    return Gmethod(name,PDd(ions),PDd(proxies),PDd(channels),
                   interferences,nblank,ndrift,ndown,
                   PAcutoff,standards,glass)
end

function PDd(arg::NamedTuple{(:P,:D,:d)})
    return PDd(arg.P,arg.D,arg.d)
end

function default_ions(name)
    m = get(_KJ["methods"],name)
    return (P=String(m.P),D=String(m.D),d=String(m.d))
end

function channels2proxies!(method::Gmethod)
    equivocal = false
    all_elements = string.(keys(_KJ["nuclides"]))
    for nuclide in (:P,:D,:d)
        channel = getproperty(method.channels,nuclide)
        matching_elements = filter(x -> occursin(x, channel), all_elements)
        if length(matching_elements) > 0
            already_found = false
            newly_found = false
            for matching_element in matching_elements
                newly_found = get_proxy_isotopes!(method,nuclide,
                                                  matching_element)
                if already_found & newly_found
                    equivocal = true
                    break
                else
                    already_found = newly_found
                end
            end
            if !already_found
                equivocal = true
            end
        else
            equivocal = true
        end
    end
    return equivocal

end
export channels2proxies!

function get_proxy_isotopes!(method::Gmethod,
                             nuclide::Symbol,
                             matching_element::AbstractString)
    channel = getproperty(method.channels,nuclide)
    all_isotopes = string.(_KJ["nuclides"][matching_element])
    matching_isotope = filter(x -> occursin(x, channel), all_isotopes)
    found = false
    if length(matching_isotope) == 1
        setproperty!(method.proxies,nuclide,
                     matching_element * matching_isotope[1])
        found = true
    end
    return found
end

function Cmethod(run::Vector{Sample},
                 standards::AbstractDict,
                 internal::Tuple;
                 nblank::Int=2)
    ch = getChannels(run)
    el = channel2element.(ch)
    elements = DataFrame(reshape(el,1,:),ch)
    return Cmethod(elements,standards,internal,nblank)
end

function getConcentrations(method::Cmethod,
                           refmat::AbstractString)
    all_concs = get(_KJ["glass"],refmat)
    elements = method.elements
    channels = names(elements)
    out = DataFrame(zeros(1, length(channels)), channels)
    for ch in channels
        out[1,ch] = all_concs[elements[1,ch]]
    end
    return out
end