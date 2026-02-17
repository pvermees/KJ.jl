function getAnchor(methodname::AbstractString,
                   refmat::AbstractString)
    standard = get(_KJ["refmat"][methodname],refmat)
    return getAnchor(standard;methodname=methodname)
end
function getAnchor(standard::IsochronRefmat;
                   methodname::AbstractString)
    t = standard.t
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
    y0 = standard.y0
    return IsochronAnchor(x0,y0,y1)
end
function getAnchor(standard::PointRefmat;kwargs...)
    return PointAnchor(standard.x,standard.y)
end
function getAnchor(standard::BiasRefmat;kwargs...)
    return BiasAnchor(standard.y0)
end
export getAnchor

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
