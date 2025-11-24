"""
Chronometer(method::String)

Initialises a geochronometer
"""
function Chronometer(method::String)
    m = get(_KJ["methods"],method)
    channels = DataFrame(par=["ion","proxy","channel"],
                         P=fill(String(m.P),3),
                         D=fill(String(m.D),3),
                         d=fill(String(m.d),3),
                         S=fill(String(m.S),3))
    nblank=2
    ndrift=2
    ndown=1
    return Chronometer(method,channels,nblank,ndrift,ndown,nothing)
end

function channelAccessor(channels::AbstractDataFrame,
                         rowname::AbstractString)
    ch = channels[findfirst(==(rowname),channels.par),2:end]
    return (P=ch.P,D=ch.D,d=ch.d,S=ch.S)
end

function getIons(chronometer::Chronometer)
    return channelAccessor(chronometer.channels,"ion")
end
export getIons

function getProxies(chronometer::Chronometer)
    return channelAccessor(chronometer.channels,"proxy")
end
export getProxies

function getChannels(chronometer::Chronometer)
    return channelAccessor(chronometer.channels,"channel")
end
export getChannels

function channelAccessor!(channels::AbstractDataFrame,
                          rowname::AbstractString;
                          P,D,d,S)
    row = findfirst(==(rowname),channels.par)
    channels[row,2:end] = [P,D,d,S]
end

function setIons!(chronometer::Chronometer;
                  P=getIons(chronometer).P,
                  D=getIons(chronometer).D,
                  d=getIons(chronometer).d,
                  S=getIons(chronometer).S)
    channelAccessor!(chronometer.channels,"ion";P,D,d,S)
end
export setIons!

function setProxies!(chronometer::Chronometer;
                     P=getIons(chronometer).P,
                     D=getIons(chronometer).D,
                     d=getIons(chronometer).d,
                     S=getIons(chronometer).S)
    channelAccessor!(chronometer.channels,"proxy";P,D,d,S)
end
export setProxies!

function setChannels!(chronometer::Chronometer;
                      P=getProxies(chronometer).P,
                      D=getProxies(chronometer).D,
                      d=getProxies(chronometer).d,
                      S=getProxies(chronometer).S)
    channelAccessor!(chronometer.channels,"channel";P,D,d,S)
end
export setChannels!
