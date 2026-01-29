function concentrations(samp::Sample,
                        method::Cmethod,
                        fit::Cfit;
                        internal::Tuple=method.internal)
    dat = swinData(samp;add_xy=true)
    sig = getSignals(dat)
    bt = polyVal(fit.blank,dat.t)
    X = getSignals(dat) .- bt
    Cs = internal[2]
    Xs = X[:,internal[1]]
    out = (X .* Cs) ./ (Xs .* fit.par)
    nms = "ppm[" .* collect(string.(values(method.elements))) .* "] from " .* names(sig)
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
    ne = length(method.elements)
    nc = 2*ne
    mat = fill(0.0,nr,nc)
    conc = nothing
    for i in eachindex(run)
        samp = run[i]
        if haskey(method.groups,samp.group)
            standard = method.groups[samp.group]
            refconcs = getConcentrations(method,standard)
            ich = method.internal[1] # ich = internal channel
            internal = (ich,refconcs[1,ich])
            conc = concentrations(samp,method,fit;
                                  internal = internal)[:,1:ne]
        else
            conc = concentrations(samp,method,fit)[:,1:ne]
        end
        mu = Statistics.mean.(eachcol(conc))
        sigma = Statistics.std.(eachcol(conc))
        nt = size(conc,1)
        mat[i,1:2:nc-1] .= mu
        mat[i,2:2:nc] .= sigma./sqrt(nt)
    end
    nms = fill("",nc)
    nms[1:2:nc-1] .= names(conc)
    nms[2:2:nc] .= "s[" .* names(conc) .* "]"
    out = hcat(DataFrame(sample=getSnames(run)),DataFrame(mat,Symbol.(nms)))
    return out
end
export concentrations
