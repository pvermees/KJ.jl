"""
internochron(run::Vector{Sample},
             channels::AbstractDict,
             blank::AbstractDataFrame,
             pars::NamedTuple;
             method::Union{AbstractString,Nothing}=nothing)

Internal isochron regression
"""
function internochron(run::Vector{Sample},
                      channels::AbstractDict,
                      blank::AbstractDataFrame,
                      pars::NamedTuple;
                      method::Union{AbstractString,Nothing}=nothing)
    ns = length(run)
    xlab = "x0"
    ylab = "y0"
    column_names = ["name", xlab, "s[" * xlab * "]", ylab, "s[" * ylab * "]", "rho"]
    out = DataFrame(hcat(fill("",ns),zeros(ns,5)),column_names)
    for i in 1:ns
        samp = run[i]
        out[i,:name] = samp.sname
        out[i,2:end] = internochron(samp,channels,blank,pars)
    end
    if isnothing(method)
        return out
    else
        return x0y02t(out,method) 
    end
end
"""
internochron(samp::Sample,
             channels::AbstractDict,
             blank::AbstractDataFrame,
             pars::NamedTuple)
"""
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
"""
internochron(Phat::AbstractVector,
             Dhat::AbstractVector,
             dhat::AbstractVector,
             vP::AbstractVector,
             vD::AbstractVector,
             vd::AbstractVector)
"""
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
    H = ForwardDiff.hessian(objective,[x0,y0])
    out = hessian2xyerr(H,[x0,y0])
    #E = covmat_internochron(x0,y0,Phat,Dhat,dhat,vP,vD,vd)
    return out
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

function SSinternochron(x0::Real,
                        y0::Real,
                        Phat::AbstractVector,
                        Dhat::AbstractVector,
                        dhat::AbstractVector,
                        vP::AbstractVector,
                        vD::AbstractVector,
                        vd::AbstractVector)
    P, D, d = internochronPDd(x0,y0,Phat,Dhat,dhat,vP,vD,vd)
    return sum(@. (P-Phat)^2/vP + (D-Dhat)^2/vd + (d-dhat)^2/vd )/2
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
    O12 = [ [ H13' H14' ]
            [ H23' H24' ] ]
    O21 = [ [ H31 H32 ]
            [ H41 H42 ] ]
    O22 = [ [ diagm(H33) diagm(H34) ]
            [ diagm(H43) diagm(H44) ] ]
    return inv( O11 - O12 * inv(O22) * O21 )
end

function x0y02t(x0y0::AbstractDataFrame,
                method::AbstractString)
    P, D, d = getPDd(method)
    xlab = "t(" * D * "/" * P * ")" 
    ylab = "(" * d * "/" * D * ")₀"
    column_names = ["name", xlab, "s[" * xlab * "]", ylab, "s[" * ylab * "]", "ρ"]
    out = DataFrame(x0y0,column_names)
    for (i,row) in enumerate(eachrow(x0y0))
        sx0y0 = row["rho"]*row["s[x0]"]*row["s[y0]"]
        E = [ [ row["s[x0]"]^2 sx0y0 ]
              [ sx0y0 row["s[y0]"]^2 ] ]
        out[i,2:end] = x0y02t(row.x0,row.y0,E,method)
    end
    return out
end
function x0y02t(x0::AbstractFloat,
                y0::AbstractFloat,
                E::Matrix,
                method::AbstractString)
    if method == "U-Pb"
        L5, L8, U58 = UPb_helper()
        init = log(1+1/x0)/L8
        objective = (par) -> york2ludwig_misfit(par,x0,y0,L5,L8,U58)
        fit = Optim.optimize(objective,[init])
        t = Optim.minimizer(fit)[1]
        dfdx0 = -y0/x0^2
        dfdy0 = 1/x0 - exp(L8*t) + 1 
        dfdt = U58*L5*exp(L5*t) - L8*y0*exp(L8*t)
        dtdx0 = -dfdx0/dfdt
        dtdy0 = -dfdy0/dfdt
    else
        lambda = _KJ["lambda"][method][1]
        t = log(1+1/x0)/lambda
        dtdx0 = -1/(lambda*x0*(1+x0))
        dtdy0 = 0.0
    end
    dy0dx0 = 0.0
    dy0dy0 = 1.0
    J = [ [ dtdx0 dtdy0 ]
          [ dy0dx0 dy0dy0 ] ]
    covmat = J * E * transpose(J)
    st = sqrt(covmat[1,1])
    sy0 = sqrt(covmat[2,2])
    rho = covmat[1,2]/(st*sy0)
    return t, st, y0, sy0, rho
end

function york2ludwig_misfit(par::AbstractVector,
                            x0::AbstractFloat,
                            y0::AbstractFloat,
                            L5::AbstractFloat,
                            L8::AbstractFloat,
                            U58::AbstractFloat)
    t = par[1]
    x = 1/(exp(L8*t)-1)
    y = U58*(exp(L5*t)-1)/(exp(L8*t)-1)
    yl = y0*(1-x/x0)
    return (yl - y)^2
end

function UPb_helper()
    L5 = _KJ["lambda"]["U235-Pb207"][1]
    L8 = _KJ["lambda"]["U238-Pb206"][1]
    fUPb = _KJ["iratio"]["U-Pb"]
    U58 = fUPb.U235/fUPb.U238
    return L5, L8, U58
end
