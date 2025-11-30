function Gmethod(name::String;
                 ions::NamedTuple{(:P,:D,:d)}=default_ions(name),
                 proxies::NamedTuple{(:P,:D,:d)}=ions,
                 channels::NamedTuple{(:P,:D,:d)}=proxies,
                 interferences::AbstractDict=Dict(),
                 nblank::Int=2,
                 ndrift::Int=2,
                 ndown::Int=1,
                 PAcutoff::Union{Nothing,Float64}=nothing,
                 anchors::AbstractDict=Dict())

    chdf = DataFrame(par=["ion","proxy","channel"],
                     P=[ions.P,proxies.P,channels.P],
                     D=[ions.D,proxies.D,channels.D],
                     d=[ions.d,proxies.d,channels.d])
    
    return Gmethod(name,chdf,interferences,nblank,ndrift,ndown,PAcutoff,anchors)
end

function Cmethod(;elements::AbstractDataFrame=DataFrame(),
                  refmats::AbstractVector=String[],
                  internal::Tuple=(nothing,nothing),
                  concentrations::AbstractDataFrame=DataFrame(),
                  nblank::Int=2)
    return Cmethod(elements,refmats,concentrations,internal,nblank)
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

function getChannels(method::Gmethod)
    return channelAccessor(method.channels,"channel")
end
function getChannels(method::Cmethod)
    return names(method.elements)
end

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
    method = method.method
    refmats = collect(keys(method.standards))
    anchors = getAnchors(method,refmats)
    method.anchors = anchors
end

function getConcentrations(method::Cmethod,
                           refmat::AbstractString)
    row = findfirst(.==(refmat),method.refmats)
    return method.concentrations[row:row,:]
end