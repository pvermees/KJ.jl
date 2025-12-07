function getAnchor(methodname::AbstractString,
                   refmat::AbstractString)
    t = get(_KJ["refmat"][methodname],refmat).type
    if t == "isochron"
        out = get_isochron_anchor(methodname,refmat)
    else
        out = get_point_anchor(methodname,refmat)
    end
    return out
end

function get_isochron_anchor(methodname::AbstractString,
                             refmat::AbstractString)
    t = get(_KJ["refmat"][methodname],refmat).tx[1]
    if methodname=="U-Pb"
        L8 = _KJ["lambda"]["U238-Pb206"][1]
        L5 = _KJ["lambda"]["U235-Pb207"][1]
        U58 = _KJ["iratio"]["U"].U235/_KJ["iratio"]["U"].U238
        x0 = 1/(exp(L8*t)-1)
        y1 = U58*(exp(L5*t)-1)/(exp(L8*t)-1)
    else
        L = _KJ["lambda"][methodname][1]
        x0 = 1/(exp(L*t)-1)
        y1 = 0.0
    end
    y0 = get(_KJ["refmat"][methodname],refmat).y0[1]
    return IsochronAnchor(x0,y0,y1)
end

function get_point_anchor(methodname::AbstractString,
                          refmat::AbstractString)
    x = get(_KJ["refmat"][methodname],refmat).tx[1]
    y = get(_KJ["refmat"][methodname],refmat).y0[1]
    return PointAnchor(x,y)
end

function getGlassAnchors(methodname::AbstractString,
                         refmats::AbstractVector)
    out = Dict()
    for refmat in refmats
        out[refmat] = get_glass_anchor(methodname,refmat)
    end
    return out
end

function get_glass_anchor(methodname::AbstractString,
                          refmat::AbstractString)
    P, D, d = getPDd(methodname)
    ratio = d * D
    y = get(_KJ["glass"],refmat)[ratio]
    return (y=y)
end
