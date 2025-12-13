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

    chdf = DataFrame(par=["ion","proxy","channel"],
                     P=[ions.P,proxies.P,channels.P],
                     D=[ions.D,proxies.D,channels.D],
                     d=[ions.d,proxies.d,channels.d])
    
    return Gmethod(name,chdf,interferences,
                   nblank,ndrift,ndown,PAcutoff,
                   standards,glass)
end

function default_ions(name)
    m = get(_KJ["methods"],name)
    return (P=String(m.P),D=String(m.D),d=String(m.d))
end

function channelAccessor(channels::AbstractDataFrame,
                         rowname::AbstractString)
    ch = channels[findfirst(==(rowname),channels.par),2:end]
    return (P=ch.P,D=ch.D,d=ch.d)
end

function getIons(method::Gmethod)
    return channelAccessor(method.channels,"ion")
end
export getIons

function getProxies(method::Gmethod)
    return channelAccessor(method.channels,"proxy")
end
export getProxies

function channelAccessor!(channels::AbstractDataFrame,
                          rowname::AbstractString;
                          P,D,d)
    row = findfirst(==(rowname),channels.par)
    channels[row,2:end] = [P,D,d]
end

function setIons!(method::Gmethod;
                  P=getIons(method).P,
                  D=getIons(method).D,
                  d=getIons(method).d)
    channelAccessor!(method.channels,"ion";P,D,d)
end
export setIons!

function setProxies!(method::Gmethod;
                     P=getIons(method).P,
                     D=getIons(method).D,
                     d=getIons(method).d)
    channelAccessor!(method.channels,"proxy";P,D,d)
end
export setProxies!

function setChannels!(method::Gmethod;
                      P=getProxies(method).P,
                      D=getProxies(method).D,
                      d=getProxies(method).d)
    channelAccessor!(method.channels,"channel";P,D,d)
end
export setChannels!

function channels2proxies!(method::Gmethod)
    equivocal = false
    all_elements = string.(keys(_KJ["nuclides"]))
    for col in eachcol(method.channels)[2:end]
        channel = col[3]
        matching_elements = filter(x -> occursin(x, channel), all_elements)
        found = false
        if length(matching_elements)==1
            found |= get_proxy_isotopes!(col,matching_elements[1])
        elseif length(matching_elements) > 1
            for matching_element in matching_elements
                found |= get_proxy_isotopes!(col,matching_element)
            end
        else
            equivocal = true
        end
        equivocal |= !found
    end
    return equivocal

end
export channels2proxies!

function get_proxy_isotopes!(col::AbstractVector,
                             matching_element::AbstractString)
    channel = col[3]
    all_isotopes = string.(_KJ["nuclides"][matching_element])
    matching_isotope = filter(x -> occursin(x, channel), all_isotopes)
    found = false
    if length(matching_isotope) > 0
        col[2] = matching_element * matching_isotope[1]
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