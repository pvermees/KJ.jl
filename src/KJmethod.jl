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
                 anchors::AbstractDict=standards2anchors(name,standards))

    chdf = DataFrame(par=["ion","proxy","channel"],
                     P=[ions.P,proxies.P,channels.P],
                     D=[ions.D,proxies.D,channels.D],
                     d=[ions.d,proxies.d,channels.d])
    
    return Gmethod(name,chdf,interferences,
                   nblank,ndrift,ndown,PAcutoff,
                   standards,anchors)
end

function standards2anchors(method::AbstractString,
                           standards::AbstractDict)
    out = Dict()
    for refmat in keys(standards)
        out[refmat] = getAnchor(method,refmat)
    end
    return out
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

function setAnchors!(method::Gmethod)
    method.anchors = Dict()
    for refmat in keys(method.standards)
        method.anchors[refmat] = getAnchor(method.name,refmat)
    end
end

function Cmethod(run::Vector{Sample},
                 standards::AbstractDict,
                 internal::Tuple;
                 nblank::Int=2)
    ch = getChannels(run)
    el = channel2element.(ch)
    elements = DataFrame(reshape(el,1,:),ch)
    concs = Dict{String,DataFrame}()
    for refmat in keys(standards)
        all_concs = get(_KJ["glass"],refmat)
        concs[refmat] = DataFrame(zeros(1, length(ch)), ch)
        for channel in names(elements)
            concs[refmat][1,channel] = all_concs[elements[1,channel]]
        end
    end
    return Cmethod(elements,standards,concs,internal,nblank)
end

function getConcentrations(method::Cmethod,
                           refmat::AbstractString)
    return method.concentrations[refmat]
end