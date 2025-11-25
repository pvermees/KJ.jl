"""
KJmethod(name::String)

Initialises a geochronometer
"""
function KJmethod(name::String)
    m = get(_KJ["methods"],name)
    channels = DataFrame(par=["ion","proxy","channel"],
                         P=fill(String(m.P),3),
                         D=fill(String(m.D),3),
                         d=fill(String(m.d),3),
                         S=fill(String(m.S),3))
    return KJmethod(name,
                    channels,
                    2,      #nblank
                    2,      #ndrift
                    1,      #ndown
                    Dict(), # standards
                    Dict()) # anchors
end

"""
setStandards!(method::KJmethod,
              standards::Dict)
"""
function setStandards!(method::KJmethod,
                       standards::Dict)
    refmats = collect(keys(standards))
    method.standards = standards
    method.anchors = getAnchors(method.name,refmats)
end
export setStandards!

function channelAccessor(channels::AbstractDataFrame,
                         rowname::AbstractString)
    ch = channels[findfirst(==(rowname),channels.par),2:end]
    return (P=ch.P,D=ch.D,d=ch.d,S=ch.S)
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
                          P,D,d,S)
    row = findfirst(==(rowname),channels.par)
    channels[row,2:end] = [P,D,d,S]
end

function setIons!(method::KJmethod;
                  P=getIons(method).P,
                  D=getIons(method).D,
                  d=getIons(method).d,
                  S=getIons(method).S)
    channelAccessor!(method.channels,"ion";P,D,d,S)
end
export setIons!

function setProxies!(method::KJmethod;
                     P=getIons(method).P,
                     D=getIons(method).D,
                     d=getIons(method).d,
                     S=getIons(method).S)
    channelAccessor!(method.channels,"proxy";P,D,d,S)
end
export setProxies!

function setChannels!(method::KJmethod;
                      P=getProxies(method).P,
                      D=getProxies(method).D,
                      d=getProxies(method).d,
                      S=getProxies(method).S)
    channelAccessor!(method.channels,"channel";P,D,d,S)
end
export setChannels!

function setAnchors!(method::KJmethod)
    method = method.method
    refmats = collect(keys(method.standards))
    anchors = getAnchors(method,refmats)
    method.anchors = anchors
end
