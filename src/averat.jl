function averat(run::Vector{Sample},
                method::Gmethod,
                fit::Gfit)
    P, D, d = getChannels(method)
    xlab = P * "/" * D
    ylab = d * "/" * D
    column_names = ["name", xlab, "s[" * xlab * "]", ylab, "s[" * ylab * "]", "rho"]
    ns = length(run)
    out = DataFrame(hcat(fill("",ns),zeros(ns,5)),column_names)
    for i in 1:ns
        samp = run[i]
        out[i,:name] = samp.sname
        out[i,2:end] = averat(samp,method,fit)
    end
    return out
end

function averat(samp::Sample,
                method::Gmethod,
                fit::Gfit)
    a = Averager(samp,method,fit)
    init = [sum(a.Phat)/sum(a.Dhat),sum(a.dhat)/sum(a.Dhat)]
    objective = (par) -> LLaverat(par[1],par[2],a)
    optimum = Optim.optimize(objective,init)
    xy = Optim.minimizer(optimum)
    H = ForwardDiff.hessian(objective,xy)
    out = hessian2xyerr(H,xy)
    return out
end
export averat

function LLaverat(x::Real,
                  y::Real,
                  a::Averager)
    Phat,Dhat,dhat,vP,vD,vd,sPD,sPd,sDd = unpack(a)
    D = @. (((dhat*vD-Dhat*sDd)*vP-Phat*sPd*vD+Dhat*sPD*sPd-dhat*sPD^2+Phat*sDd*sPD)*y+((Phat*vD-Dhat*sPD)*vd-dhat*sPd*vD+Dhat*sDd*sPd+dhat*sDd*sPD-Phat*sDd^2)*x+(Dhat*vP-Phat*sPD)*vd-dhat*sDd*vP-Dhat*sPd^2+(dhat*sPD+Phat*sDd)*sPd)/((vD*vP-sPD^2)*y^2+((2*sDd*sPD-2*sPd*vD)*x-2*sDd*vP+2*sPD*sPd)*y+(vD*vd-sDd^2)*x^2+(2*sDd*sPd-2*sPD*vd)*x+vP*vd-sPd^2)
    SS = @. (D*y-dhat)*(((vD*vP-sPD^2)*(D*y-dhat))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((sDd*sPD-sPd*vD)*(D*x-Phat))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((D-Dhat)*(sPD*sPd-sDd*vP))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD)))+(D-Dhat)*(((sPD*sPd-sDd*vP)*(D*y-dhat))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((sDd*sPd-sPD*vd)*(D*x-Phat))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((D-Dhat)*(vP*vd-sPd^2))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD)))+(D*x-Phat)*(((sDd*sPD-sPd*vD)*(D*y-dhat))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((vD*vd-sDd^2)*(D*x-Phat))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((D-Dhat)*(sDd*sPd-sPD*vd))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD)))
    return sum(@. SS/2)
end

function Averager(samp::Sample,
                  method::Gmethod,
                  fit::Gfit)
    a = atomic(samp,method,fit)
    mat = hcat(a.P,a.D,a.d)
    E = df2cov(mat)
    return Averager(a.P,a.D,a.d,E[1,1],E[2,2],E[3,3],E[1,2],E[1,3],E[2,3])
end