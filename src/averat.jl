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
        out[i,2:end] = averat(samp,channels,blank,pars)
    end
    return out
end
function averat(samp::Sample,
                channels::AbstractDict,
                blank::AbstractDataFrame,
                pars::NamedTuple)
    Phat, Dhat, dhat = atomic(samp,channels,blank,pars)
    return averat(Phat,Dhat,dhat)
end
function averat(Phat::AbstractVector,
                Dhat::AbstractVector,
                dhat::AbstractVector)
    O = O_timeseries(Phat,Dhat,dhat)
    init = [sum(Phat)/sum(Dhat),sum(dhat)/sum(Dhat)]
    objective = (par) -> LLaverat(par[1],par[2],
                                  Phat,Dhat,dhat,O)
    fit = Optim.optimize(objective,init)
    x, y = Optim.minimizer(fit)
    H = ForwardDiff.hessian(objective,[x,y])
    out = hessian2xyerr(H,[x,y])
    return out
end
export averat

function LLaverat(x::Real,
                  y::Real,
                  Phat::AbstractVector,
                  Dhat::AbstractVector,
                  dhat::AbstractVector,
                  O::Matrix)
    D = averatD(x,y,Phat,Dhat,dhat,O)
    return sum(@. ((D*y-dhat)*(O[3,3]*(D*y-dhat)+O[3,1]*(D*x-Phat)+(D-Dhat)*O[3,2])+(D-Dhat)*(O[2,3]*(D*y-dhat)+O[2,1]*(D*x-Phat)+(D-Dhat)*O[2,2])+(D*x-Phat)*(O[1,3]*(D*y-dhat)+O[1,1]*(D*x-Phat)+(D-Dhat)*O[1,2]))/2 )
end

function averatD(x::Real,
                 y::Real,
                 Phat::AbstractVector,
                 Dhat::AbstractVector,
                 dhat::AbstractVector,
                 O::Matrix)
    return @. ((2*O[3,3]*dhat+(O[3,1]+O[1,3])*Phat+Dhat*O[3,2]+Dhat*O[2,3])*y+((O[3,1]+O[1,3])*dhat+2*O[1,1]*Phat+Dhat*O[2,1]+Dhat*O[1,2])*x+(O[3,2]+O[2,3])*dhat+(O[2,1]+O[1,2])*Phat+2*Dhat*O[2,2])/(2*O[3,3]*y^2+((2*O[3,1]+2*O[1,3])*x+2*O[3,2]+2*O[2,3])*y+2*O[1,1]*x^2+(2*O[2,1]+2*O[1,2])*x+2*O[2,2])
end
