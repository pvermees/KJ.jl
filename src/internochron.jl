function internochron(run::Vector{Sample},
                      method::AbstractString,
                      channels::AbstractDict,
                      blank::AbstractDataFrame,
                      pars::NamedTuple)
    ns = length(run)
    P, D, d = getPDd(method)
    xlab = P * "/" * D
    ylab = d * "/" * D
    column_names = ["name", xlab, "s[" * xlab * "]", ylab, "s[" * ylab * "]", "rho"]
    out = DataFrame(hcat(fill("",ns),zeros(ns,5)),column_names)
    for i in 1:ns
        samp = run[i]
        out[i,:name] = samp.sname
        out[i,2:end] = internochron(samp,channels,blank,pars)
    end
    return out
end
function internochron(samp::Sample,
                      channels::AbstractDict,
                      blank::AbstractDataFrame,
                      pars::NamedTuple)
    Phat, Dhat, dhat = atomic(samp,channels,blank,pars)
    vP = var_timeseries(Phat)
    vD = var_timeseries(Dhat)
    vd = var_timeseries(dhat)
    return internochron(Phat,Dhat,dhat,vP,vD,vd)
end
function internochron(Phat::AbstractVector,
                      Dhat::AbstractVector,
                      dhat::AbstractVector,
                      vP::AbstractVector,
                      vD::AbstractVector,
                      vd::AbstractVector)
    init = init_internochron(Phat,Dhat,dhat)
    objective = (par) -> SSinternochron(par[1],par[2],
                                        Phat,Dhat,dhat,
                                        vP,vD,vd)
    fit = Optim.optimize(objective,init)
    x0, y0 = Optim.minimizer(fit)
    E = covmat_internochron(x0,y0,Phat,Dhat,dhat,vP,vD,vd)
    sx0 = E[1,1]>0 ? sqrt(E[1,1]) : NaN
    sy0 = E[2,2]>0 ? sqrt(E[2,2]) : NaN
    rx0y0 = E[1,2]/(sx0*sy0)
    return [x0 sx0 y0 sy0 rx0y0]
end
export internochron

function init_internochron(Phat::AbstractVector,
                           Dhat::AbstractVector,
                           dhat::AbstractVector)
    minx = minimum(Phat./Dhat)
    maxx = maximum(Phat./Dhat)
    miny = minimum(dhat./Dhat)
    maxy = maximum(dhat./Dhat)
    slope = (maxy-miny)/(maxx-minx)
    y0i = maxy+slope*minx
    x0i = y0i/slope
    return [x0i,y0i]
end

function SSinternochron(x0::AbstractFloat,
                        y0::AbstractFloat,
                        Phat::AbstractVector,
                        Dhat::AbstractVector,
                        dhat::AbstractVector,
                        vP::AbstractVector,
                        vD::AbstractVector,
                        vd::AbstractVector)
    P, D, d = internochronPDd(x0,y0,Phat,Dhat,dhat,vP,vD,vd)
    return sum(@. (P-Phat)^2/vP + (D-Dhat)^2/vd + (d-dhat)^2/vd )
end

function internochronPd(x0,y0,Phat,Dhat,dhat,vP,vD,vd)
    P = @. ((Phat*vD*x0^2+Dhat*vP*x0)*y0^2-dhat*vP*x0*y0+Phat*vd*x0^2)/((vD*x0^2+vP)*y0^2+vd*x0^2)
    d = @. ((dhat*vD*x0^2+dhat*vP)*y0^2+(Dhat*vd*x0^2-Phat*vd*x0)*y0)/((vD*x0^2+vP)*y0^2+vd*x0^2)
    return P, d
end

function internochronPDd(x0,y0,Phat,Dhat,dhat,vP,vD,vd)
    P, d = internochronPd(x0,y0,Phat,Dhat,dhat,vP,vD,vd)
    D = @. d/y0 + P/x0
    return P, D, d
end

function covmat_internochron(x0,y0,Phat,Dhat,dhat,vP,vD,vd)
    P, D, d = internochronPDd(x0,y0,Phat,Dhat,dhat,vP,vD,vd)
    H11 = sum(@. P^2/(vD*x0^4)-(2*P*((-d/y0)-P/x0+Dhat))/(vD*x0^3) )
    H12 = sum(@. (P*d)/(vD*x0^2*y0^2) )
    H13 = @. ((-d/y0)-P/x0+Dhat)/(vD*x0^2)-P/(vD*x0^3)
    H14 = @. -P/(vD*x0^2*y0)
    H21 = sum(@. (P*d)/(vD*x0^2*y0^2) )
    H22 = sum(@. d^2/(vD*y0^4)-(2*d*((-d/y0)-P/x0+Dhat))/(vD*y0^3) )
    H23 = @. -d/(vD*x0*y0^2)
    H24 = @. ((-d/y0)-P/x0+Dhat)/(vD*y0^2)-d/(vD*y0^3)
    H31 = @. ((2*((-d/y0)-P/x0+Dhat))/(vD*x0^2)-(2*P)/(vD*x0^3))/2
    H32 = @. -d/(vD*x0*y0^2)
    H33 = @. (2/(vD*x0^2)+2/vP)/2
    H34 = @. 1/(vD*x0*y0)
    H41 = @. -P/(vD*x0^2*y0)
    H42 = @. ((2*((-d/y0)-P/x0+Dhat))/(vD*y0^2)-(2*d)/(vD*y0^3))/2
    H43 = @. 1/(vD*x0*y0)
    H44 = @. (2/(vD*y0^2)+2/vd)/2
    O11 = [ [H11 H12]
            [H21 H22] ]
    O12 = [ [ H13 H14 ]
            [ H23 H24 ] ]'
    O21 = [ [ H31 H32 ]
            [ H41 H42 ] ]
    O22 = [ [ diagm(H33) diagm(H34) ]
            [ diagm(H43) diagm(H44) ] ]
    return inv( O11 - O12 * inv(O22) * O21 )
end
