function Gmethod(name::String;
                 ions::NamedTuple{(:P,:D,:d)}=default_ions(name),
                 channels::NamedTuple{(:P,:D,:d)}=ions,
                 proxies::NamedTuple{(:P,:D,:d)}=channels2proxies(channels),
                 refmats::AbstractDict=Dict{String,String}(),
                 standards::AbstractVector=collect(keys(refmats)),
                 bias::AbstractDict=Dict{String,Vector{String}}(),
                 fractionation::Fractionation=Fractionation(ions,proxies,channels,standards,bias),
                 interference::Interference=Interference(),
                 nblank::Int=2,
                 ndrift::Int=2,
                 ndown::Int=1,
                 nbias::Int=1,
                 PAcutoff::Union{Nothing,Float64}=nothing)
    return Gmethod(name,refmats,fractionation,interference,
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
    if length(matching_elements) > 0
        already_found = false
        for matching_element in matching_elements
            proxy = get_proxy_isotope(channel,matching_element)
            if !isnothing(proxy)
                if already_found
                    return nothing
                else
                    already_found = true
                end
            end
        end
    end
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
                 refmats::AbstractDict=Dict{String,String}(),
                 internal::Tuple=(nothing,nothing),
                 nblank::Int=2)
    ch = getChannels(run)
    el = channel2element.(ch)
    elements = NamedTuple{Tuple(Symbol.(ch))}(Tuple(el))
    return Cmethod(elements,refmats,internal,nblank)
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