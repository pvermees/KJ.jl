"""
averat

Average the 'atomic' isotopic ratios for a sample

# Returns

- a dataframe of P/D and d/D-ratios with their standard errors and error correlations

# Arguments

See [`process!`](@ref).
"""
function averat(run::Vector{Sample},
                channels::AbstractDict,
                blank::AbstractDataFrame,
                pars::NamedTuple;
                dt::Union{AbstractDict,Nothing}=nothing,
                dead::AbstractFloat=0.0,
                method=nothing)
    ns = length(run)
    if isnothing(method)
        xlab = "x"
        ylab = "y"
    else
        P, D, d = getPDd(method)
        xlab = P * "/" * D
        ylab = d * "/" * D
    end
    column_names = ["name", xlab, "s[" * xlab * "]", ylab, "s[" * ylab * "]", "rho"]
    out = DataFrame(hcat(fill("",ns),zeros(ns,5)),column_names)
    for i in 1:ns
        samp = run[i]
        out[i,:name] = samp.sname
        out[i,2:end] = averat(samp,channels,blank,pars;
                              dt=dt,dead=dead)
    end
    return out
end
function averat(samp::Sample,
                channels::AbstractDict,
                blank::AbstractDataFrame,
                pars::NamedTuple;
                dt::Union{AbstractDict,Nothing}=nothing,
                dead::AbstractFloat=0.0)
    Phat, Dhat, dhat = atomic(samp,channels,blank,pars;
                              dt=dt,dead=dead)
    return averat(Phat,Dhat,dhat;dt=dt,dead=dead)
end
function averat(Phat::AbstractVector,
                Dhat::AbstractVector,
                dhat::AbstractVector;
                dt::Union{AbstractDict,Nothing}=nothing,
                dead::AbstractFloat=0.0)
    if isnothing(dt)
        vP = var_timeseries(Phat)
        vD = var_timeseries(Dhat)
        vd = var_timeseries(dhat)
    else
        vP = var_cps(Phat,dt,dead)
        vD = var_cps(Dhat,dt,dead)
        vd = var_cps(dhat,dt,dead)
    end
    init = [sum(Phat)/sum(Dhat),sum(dhat)/sum(Dhat)]
    objective = (par) -> SSaverat(par[1],par[2],
                                  Phat,Dhat,dhat,
                                  vP,vD,vd)
    fit = Optim.optimize(objective,init)
    x, y = Optim.minimizer(fit)
    E = covmat_averat(x,y,Phat,Dhat,dhat,vP,vD,vd)
    sx = sqrt(E[1,1])
    sy = sqrt(E[2,2])
    rxy = E[1,2]/(sx*sy)
    return [x sx y sy rxy]
end
export averat

function SSaverat(x::AbstractFloat,
                  y::AbstractFloat,
                  Phat::AbstractVector,
                  Dhat::AbstractVector,
                  dhat::AbstractVector,
                  vP::AbstractVector,
                  vD::AbstractVector,
                  vd::AbstractVector)
    S = averatS(x,y,Phat,Dhat,dhat,vP,vD,vd)
    dP = @. Phat - S*x/(1+x+y)
    dd = @. dhat - S*y/(1+x+y)
    dD = @. Dhat - S/(1+x+y)
    return sum(@. dP^2/vP + dD^2/vD + dd^2/vd )
end

# block matrix inversion of the Hessian matrix
function covmat_averat(x::AbstractFloat,
                       y::AbstractFloat,
                       Phat::AbstractVector,
                       Dhat::AbstractVector,
                       dhat::AbstractVector,
                       vP::AbstractVector,
                       vD::AbstractVector,
                       vd::AbstractVector)
    S = averatS(x,y,Phat,Dhat,dhat,vP,vD,vd)
    O11 = sum(@. ((2*(S/(y+x+1)-(S*x)/(y+x+1)^2)^2)/vP+(4*S*y*((S*y)/(y+x+1)-dhat))/(vd*(y+x+1)^3)+(2*((2*S*x)/(y+x+1)^3-(2*S)/(y+x+1)^2)*((S*x)/(y+x+1)-Phat))/vP+(4*S*(S/(y+x+1)-Dhat))/(vD*(y+x+1)^3)+(2*S^2*y^2)/(vd*(y+x+1)^4)+(2*S^2)/(vD*(y+x+1)^4))/2 )
    O12 = sum(@. ((-(2*S*((S*y)/(y+x+1)-dhat))/(vd*(y+x+1)^2))+(4*S*y*((S*y)/(y+x+1)-dhat))/(vd*(y+x+1)^3)+(2*((2*S*x)/(y+x+1)^3-S/(y+x+1)^2)*((S*x)/(y+x+1)-Phat))/vP-(2*S*y*(S/(y+x+1)-(S*y)/(y+x+1)^2))/(vd*(y+x+1)^2)-(2*S*x*(S/(y+x+1)-(S*x)/(y+x+1)^2))/(vP*(y+x+1)^2)+(4*S*(S/(y+x+1)-Dhat))/(vD*(y+x+1)^3)+(2*S^2)/(vD*(y+x+1)^4))/2 )
    O21 = sum(@. ((2*((2*S*y)/(y+x+1)^3-S/(y+x+1)^2)*((S*y)/(y+x+1)-dhat))/vd-(2*S*((S*x)/(y+x+1)-Phat))/(vP*(y+x+1)^2)+(4*S*x*((S*x)/(y+x+1)-Phat))/(vP*(y+x+1)^3)-(2*S*y*(S/(y+x+1)-(S*y)/(y+x+1)^2))/(vd*(y+x+1)^2)-(2*S*x*(S/(y+x+1)-(S*x)/(y+x+1)^2))/(vP*(y+x+1)^2)+(4*S*(S/(y+x+1)-Dhat))/(vD*(y+x+1)^3)+(2*S^2)/(vD*(y+x+1)^4))/2 )
    O22 = sum(@. ((2*(S/(y+x+1)-(S*y)/(y+x+1)^2)^2)/vd+(2*((2*S*y)/(y+x+1)^3-(2*S)/(y+x+1)^2)*((S*y)/(y+x+1)-dhat))/vd+(4*S*x*((S*x)/(y+x+1)-Phat))/(vP*(y+x+1)^3)+(4*S*(S/(y+x+1)-Dhat))/(vD*(y+x+1)^3)+(2*S^2*x^2)/(vP*(y+x+1)^4)+(2*S^2)/(vD*(y+x+1)^4))/2)
    O13 = @. ((-(2*y*((S*y)/(y+x+1)-dhat))/(vd*(y+x+1)^2))+(2*(1/(y+x+1)-x/(y+x+1)^2)*((S*x)/(y+x+1)-Phat))/vP+(2*x*(S/(y+x+1)-(S*x)/(y+x+1)^2))/(vP*(y+x+1))-(2*(S/(y+x+1)-Dhat))/(vD*(y+x+1)^2)-(2*S*y^2)/(vd*(y+x+1)^3)-(2*S)/(vD*(y+x+1)^3))/2
    O23 = @. ((2*(1/(y+x+1)-y/(y+x+1)^2)*((S*y)/(y+x+1)-dhat))/vd-(2*x*((S*x)/(y+x+1)-Phat))/(vP*(y+x+1)^2)+(2*y*(S/(y+x+1)-(S*y)/(y+x+1)^2))/(vd*(y+x+1))-(2*(S/(y+x+1)-Dhat))/(vD*(y+x+1)^2)-(2*S*x^2)/(vP*(y+x+1)^3)-(2*S)/(vD*(y+x+1)^3))/2
    O31 = @. ((-(2*y*((S*y)/(y+x+1)-dhat))/(vd*(y+x+1)^2))+(2*((S*x)/(y+x+1)-Phat))/(vP*(y+x+1))-(2*x*((S*x)/(y+x+1)-Phat))/(vP*(y+x+1)^2)+(2*x*(S/(y+x+1)-(S*x)/(y+x+1)^2))/(vP*(y+x+1))-(2*(S/(y+x+1)-Dhat))/(vD*(y+x+1)^2)-(2*S*y^2)/(vd*(y+x+1)^3)-(2*S)/(vD*(y+x+1)^3))/2
    O32 = @. ((2*((S*y)/(y+x+1)-dhat))/(vd*(y+x+1))-(2*y*((S*y)/(y+x+1)-dhat))/(vd*(y+x+1)^2)-(2*x*((S*x)/(y+x+1)-Phat))/(vP*(y+x+1)^2)+(2*y*(S/(y+x+1)-(S*y)/(y+x+1)^2))/(vd*(y+x+1))-(2*(S/(y+x+1)-Dhat))/(vD*(y+x+1)^2)-(2*S*x^2)/(vP*(y+x+1)^3)-(2*S)/(vD*(y+x+1)^3))/2
    O33 = @. ((2*y^2)/(vd*(y+x+1)^2)+(2*x^2)/(vP*(y+x+1)^2)+2/(vD*(y+x+1)^2))/2
    H11 = [ [O11 O12]
            [O21 O22] ]
    H12 = [ O13' ; O23' ]
    H21 = [ O31 O32 ]
    H22 = diagm( O33 )
    return inv( H11 - H12 * inv(H22) * H21 )
end

function averatS(x::AbstractFloat,
                 y::AbstractFloat,
                 Phat::AbstractVector,
                 Dhat::AbstractVector,
                 dhat::AbstractVector,
                 vP::AbstractVector,
                 vD::AbstractVector,
                 vd::AbstractVector)
    return @. (dhat*vD*vP*y^2+((Phat*vD*vd+dhat*vD*vP)*x+Dhat*vP*vd+dhat*vD*vP)*y +Phat*vD*vd*x^2+(Dhat*vP+Phat*vD)*vd*x+Dhat*vP*vd)/(vD*vP*y^2+vD*vd*x^2+vP*vd)
end
