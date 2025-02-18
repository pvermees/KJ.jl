function internochron(run::Vector{Sample},
                      channels::AbstractDict,
                      blank::AbstractDataFrame,
                      pars::NamedTuple;
                      method::Union{AbstractString,Nothing}=nothing)
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
        out[i,2:end] = internochron(samp,channels,blank,pars)
    end
    return out
end
function internochron(samp::Sample,
                      channels::AbstractDict,
                      blank::AbstractDataFrame,
                      pars::NamedTuple)
    Phat, Dhat, dhat = atomic(samp,channels,blank,pars)
    vPhat = var_timeseries(Phat)
    vDhat = var_timeseries(Dhat)
    vdhat = var_timeseries(dhat)
    return internochron(Phat,Dhat,dhat,vPhat,vDhat,vdhat)
end
function internochron(Phat::AbstractVector,
                      Dhat::AbstractVector,
                      dhat::AbstractVector,
                      vPhat::AbstractVector,
                      vDhat::AbstractVector,
                      vdhat::AbstractVector)
    init = init_internochron(Phat,Dhat,dhat)
    objective = (par) -> SSinternochron(par[1],par[2],
                                        Phat,Dhat,dhat,
                                        vPhat,vDhat,vdhat)
    fit = Optim.optimize(objective,init)
    x0hat, y0hat = Optim.minimizer(fit)
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
                        vPhat::AbstractVector,
                        vDhat::AbstractVector,
                        vdhat::AbstractVector)
    P, D, d = internochronPDd(x0,y0,Phat,Dhat,dhat,vPhat,vDhat,vdhat)
    return sum(@. (P-Phat)^2/vPhat + (D-Dhat)^2/vDhat + (d-dhat)^2/vdhat )
end

function internochronPd(x0,y0,Phat,Dhat,dhat,vPhat,vDhat,vdhat)
    P = @. ((Phat*vDhat*x0^2+Dhat*vPhat*x0)*y0^2-dhat*vPhat*x0*y0+Phat*vdhat*x0^2)/((vDhat*x0^2+vPhat)*y0^2+vdhat*x0^2)
    d = @. ((dhat*vDhat*x0^2+dhat*vPhat)*y0^2+(Dhat*vdhat*x0^2-Phat*vdhat*x0)*y0)/((vDhat*x0^2+vPhat)*y0^2+vdhat*x0^2)
    return P, d
end

function internochronPDd(x0,y0,Phat,Dhat,dhat,vPhat,vDhat,vdhat)
    P, d = internochronPd(x0,y0,Phat,Dhat,dhat,vPhat,vDhat,vdhat)
    D = @. d/y0 + P/x0
    return P, D, d
end
