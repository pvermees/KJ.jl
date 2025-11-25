"""
KJmethod(name::String)

Initialises a geochronometer
"""
function KJmethod(name::String)
    m = get(_KJ["methods"],name)
    channels = DataFrame(par=["ion","proxy","channel"],
                         P=fill(String(m.P),3),
                         D=fill(String(m.D),3),
                         d=fill(String(m.d),3))
    return KJmethod(name,
                    channels,
                    Dict(),  # interferences
                    2,       # nblank
                    2,       # ndrift
                    1,       # ndown
                    nothing, # PAcutoff
                    Dict())  # anchors
end

function channelAccessor(channels::AbstractDataFrame,
                         rowname::AbstractString)
    ch = channels[findfirst(==(rowname),channels.par),2:end]
    return (P=ch.P,D=ch.D,d=ch.d)
end

function getIons(method::KJmethod)
    return channelAccessor(method.channels,"ion")
end
export getIons

function getProxies(method::KJmethod)
    return channelAccessor(method.channels,"proxy")
end
export getProxies

function getChannels(method::KJmethod)
    return channelAccessor(method.channels,"channel")
end
export getChannels

function channelAccessor!(channels::AbstractDataFrame,
                          rowname::AbstractString;
                          P,D,d)
    row = findfirst(==(rowname),channels.par)
    channels[row,2:end] = [P,D,d]
end

function setIons!(method::KJmethod;
                  P=getIons(method).P,
                  D=getIons(method).D,
                  d=getIons(method).d)
    channelAccessor!(method.channels,"ion";P,D,d)
end
export setIons!

function setProxies!(method::KJmethod;
                     P=getIons(method).P,
                     D=getIons(method).D,
                     d=getIons(method).d)
    channelAccessor!(method.channels,"proxy";P,D,d)
end
export setProxies!

function setChannels!(method::KJmethod;
                      P=getProxies(method).P,
                      D=getProxies(method).D,
                      d=getProxies(method).d)
    channelAccessor!(method.channels,"channel";P,D,d)
end
export setChannels!

function setAnchors!(method::KJmethod)
    method = method.method
    refmats = collect(keys(method.standards))
    anchors = getAnchors(method,refmats)
    method.anchors = anchors
end
