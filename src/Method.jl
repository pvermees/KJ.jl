"""
Method(method::AbstractString)

Initialises a
"""
function Method(method::AbstractString)
    m = _KJ["methods"][method,:]
    ions = (P=m.P,D=m.D,d=m.D,S=m.S)
    return PDdS(method,ions,nothing,nothing)
end

function setProxies!(method::Method;
                     P::AbstractString=method.ions.P,
                     D::AbstractString=method.ions.D,
                     d::AbstractString=method.ions.d
                     S::AbstractString=method.ions.S)
    proxies = (P=m.P,D=m.D,d=m.D,S=m.S)
    method.proxies = proxies
end

function setChannels!(method::Method;
                      P::AbstractString=ifelse(isnothing(method.proxies),method.ions.P,method.proxies.P),
                      D::AbstractString=ifelse(isnothing(method.proxies),method.ions.D,method.proxies.D),
                      d::AbstractString=ifelse(isnothing(method.proxies),method.ions.d,method.proxies.d),
                      S::AbstractString=ifelse(isnothing(method.proxies),method.ions.d,method.proxies.S))
    proxies = (P=m.P,D=m.D,d=m.D,S=m.S)
    method.proxies = proxies
end
