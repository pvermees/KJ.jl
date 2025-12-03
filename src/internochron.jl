function internochron(run::Vector{Sample},
                      method::Gmethod,
                      fit::Gfit)
    ns = length(run)
    xlab = "x0"
    ylab = "y0"
    column_names = ["name", xlab, "s[" * xlab * "]", ylab, "s[" * ylab * "]", "rho"]
    out = DataFrame(hcat(fill("",ns),zeros(ns,5)),column_names)
    for i in 1:ns
        samp = run[i]
        out[i,:name] = samp.sname
        out[i,2:end] = internochron(samp,method,fit)
    end
    return x0y02t(out,method)
end

function internochron(samp::Sample,
                      method::Gmethod,
                      fit::Gfit)
    a = Averager(samp,method,fit)
    init = init_internochron(a)
    objective = (par) -> LLinternochron(par[1],par[2],a)
    fit = Optim.optimize(objective,init)
    x0, y0 = Optim.minimizer(fit)
    H = ForwardDiff.hessian(objective,[x0,y0])
    out = hessian2xyerr(H,[x0,y0])
    return out
end
export internochron

function init_internochron(a::Averager)
    minx = minimum(a.Phat./a.Dhat)
    maxx = maximum(a.Phat./a.Dhat)
    miny = minimum(a.dhat./a.Dhat)
    maxy = maximum(a.dhat./a.Dhat)
    slope = (maxy-miny)/(maxx-minx)
    y0i = maxy+slope*minx
    x0i = y0i/slope
    return [x0i,y0i]
end

function LLinternochron(x0::Real,
                        y0::Real,
                        a::Averager)
    Phat,Dhat,dhat,vP,vD,vd,sPD,sPd,sDd = unpack(a)
    P = @. (((Phat*vD-Dhat*sPD)*x0^2+(Dhat*vP-Phat*sPD)*x0)*y0^2+((Dhat*sPd+dhat*sPD-2*Phat*sDd)*x0^2+(Phat*sPd-dhat*vP)*x0)*y0+(Phat*vd-dhat*sPd)*x0^2)/((vD*x0^2-2*sPD*x0+vP)*y0^2+(2*sPd*x0-2*sDd*x0^2)*y0+vd*x0^2)
    d = @. (((dhat*vD-Dhat*sDd)*x0^2+(Dhat*sPd-2*dhat*sPD+Phat*sDd)*x0+dhat*vP-Phat*sPd)*y0^2+((Dhat*vd-dhat*sDd)*x0^2+(dhat*sPd-Phat*vd)*x0)*y0)/((vD*x0^2-2*sPD*x0+vP)*y0^2+(2*sPd*x0-2*sDd*x0^2)*y0+vd*x0^2)
    SS = @. (-(d/y0)-P/x0+Dhat)*(((vP*vd-sPd^2)*(-(d/y0)-P/x0+Dhat))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((Phat-P)*(sDd*sPd-sPD*vd))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((dhat-d)*(sPD*sPd-sDd*vP))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD)))+(Phat-P)*(((sDd*sPd-sPD*vd)*(-(d/y0)-P/x0+Dhat))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((Phat-P)*(vD*vd-sDd^2))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((dhat-d)*(sDd*sPD-sPd*vD))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD)))+(dhat-d)*(((sPD*sPd-sDd*vP)*(-(d/y0)-P/x0+Dhat))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((dhat-d)*(vD*vP-sPD^2))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((Phat-P)*(sDd*sPD-sPd*vD))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD)))
    return sum(@. SS/2)
end

function x0y02t(x0y0::AbstractDataFrame,
                method::Gmethod)
    P, D, d = getChannelsDict(method)
    xlab = "t(" * D * "/" * P * ")" 
    ylab = "(" * d * "/" * D * ")₀"
    column_names = ["name", xlab, "s[" * xlab * "]", ylab, "s[" * ylab * "]", "ρ"]
    out = DataFrame(x0y0,column_names)
    for (i,row) in enumerate(eachrow(x0y0))
        sx0y0 = row["rho"]*row["s[x0]"]*row["s[y0]"]
        E = [ [ row["s[x0]"]^2 sx0y0 ]
              [ sx0y0 row["s[y0]"]^2 ] ]
        ty0 = x0y02t(row.x0,row.y0,E,method.name)
        out[i,2:end] = values(ty0)
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
    return (t=t, st=st, y0=y0, sy0=sy0, rho=rho)
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
    U58 = _KJ["iratio"]["U"].U235/_KJ["iratio"]["U"].U238
    return L5, L8, U58
end
