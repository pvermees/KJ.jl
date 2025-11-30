function concentrations(samp::Sample,
                        method::Cmethod,
                        fit::Cfit)
    dat = swinData(samp;add_xy=true)
    sig = getSignals(dat)
    bt = polyVal(fit.blank,dat.t)
    X = getSignals(dat) .- bt
    Cs = method.internal[2]
    Xs = X[:,method.internal[1]]
    out = (X .* Cs) ./ (Xs .* fit.par)
    nms = "ppm[" .* Vector(method.elements[1,:]) .* "] from " .* names(sig)
    if "x" in names(dat) && "y" in names(dat)
        out.x = dat.x
        out.y = dat.y
        append!(nms,["x","y"])
    end
    rename!(out,Symbol.(nms))
    return out
end
function concentrations(run::Vector{Sample},
                        method::Cmethod,
                        fit::Cfit)
    nr = length(run)
    nc = 2*ncol(method.elements)
    mat = fill(0.0,nr,nc)
    conc = nothing
    for i in eachindex(run)
        samp = run[i]
        conc = concentrations(samp,method,fit)
        mu = Statistics.mean.(eachcol(conc))
        sigma = Statistics.std.(eachcol(conc))
        mat[i,1:2:nc-1] .= mu
        mat[i,2:2:nc] .= sigma
    end
    nms = fill("",nc)
    nms[1:2:nc-1] .= names(conc)
    nms[2:2:nc] .= "s[" .* names(conc) .* "]"
    out = hcat(DataFrame(sample=getSnames(run)),DataFrame(mat,Symbol.(nms)))
    return out
end
export concentrations
