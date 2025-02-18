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
    D = avaratD(x,y,Phat,Dhat,dhat,vP,vD,vd)
    sum(@. (D*y-dhat)^2/vd+(D*x-Phat)^2/vP+(D-Dhat)^2/vD )
end

# block matrix inversion of the Hessian matrix
function covmat_averat(x,y,Phat,Dhat,dhat,vP,vD,vd)
    D = averatD(x,y,Phat,Dhat,dhat,vP,vD,vd)
    H11 = [ [sum(@. D^2/vP) 0]
            [0 sum(@. D^2/vd)] ]
    H12 = [ (@. (D*x-Phat)/vP+(D*x)/vP)' ;
            (@. (D*y-dhat)/vd+(D*y)/vd)' ]
    H21 = [ (@. ((2*(D*x-Phat))/vP+(2*D*x)/vP)/2)' ;
            (@. ((2*(D*y-dhat))/vd+(2*D*y)/vd)/2)' ]'
    H22 = diagm( @. ((2*y^2)/vd+(2*x^2)/vP+2/vD)/2 )
    return inv( H11 - H12 * inv(H22) * H21 )
end

function averatD(x,y,Phat,Dhat,dhat,vP,vD,vd)
    return @. (dhat*vD*vP*y+Phat*vD*vd*x+Dhat*vP*vd)/(vD*vP*y^2+vD*vd*x^2+vP*vd)
end
