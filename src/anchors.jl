function getAnchors(method::AbstractString,
                    refmats::AbstractVector)
    out = Dict()
    for refmat in refmats
        out[refmat] = getAnchor(method,refmat)
    end
    return out
end

function getAnchor(method::AbstractString,
                   refmat::AbstractString)
    t = get(_KJ["refmat"][method],refmat).type
    if t == "isochron"
        out = get_isochron_anchor(method,refmat)
    else
        out = get_point_anchor(method,refmat)
    end
    return out
end

function get_isochron_anchor(method::AbstractString,
                             refmat::AbstractString)
    t = get(_KJ["refmat"][method],refmat).tx[1]
    if method=="U-Pb"
        L8 = _KJ["lambda"]["U238-Pb206"][1]
        L5 = _KJ["lambda"]["U235-Pb207"][1]
        U58 = _KJ["iratio"]["U"].U235/_KJ["iratio"]["U"].U238
        x0 = 1/(exp(L8*t)-1)
        y1 = U58*(exp(L5*t)-1)/(exp(L8*t)-1)
    else
        L = _KJ["lambda"][method][1]
        x0 = 1/(exp(L*t)-1)
        y1 = 0.0
    end
    y0 = get(_KJ["refmat"][method],refmat).y0[1]
    return (x0=x0,y0=y0,y1=y1)
end

function is_isochron_anchor(anchor::NamedTuple)
    return all(in(keys(anchor)), [:x0,:y0,:y1])
end

function get_point_anchor(method::AbstractString,
                          refmat::AbstractString)
    x = get(_KJ["refmat"][method],refmat).tx[1]
    y = get(_KJ["refmat"][method],refmat).y0[1]
    return (x=x,y=y)
end

function is_point_anchor(anchor::NamedTuple)
    k = keys(anchor)
    return all(in(keys(anchor)), [:x,:y])
end

function getGlassAnchors(method::AbstractString,
                         refmats::AbstractVector)
    out = Dict()
    for refmat in refmats
        out[refmat] = get_glass_anchor(method,refmat)
    end
    return out
end

function get_glass_anchor(method::AbstractString,
                          refmat::AbstractString)
    P, D, d = getPDd(method)
    ratio = d * D
    y = get(_KJ["glass"],refmat)[ratio]
    return (y=y)
end
