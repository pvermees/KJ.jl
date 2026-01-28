function fractionation!(fit::Gfit,
                        method::Gmethod,
                        run::Vector{Sample};
                        verbose::Bool=false)

    # extract the grouped data for the SS function from the run
    cruncher_groups = Dict()
    for group in method.standards
        standard = method.groups[group]
        anchor = getAnchor(method.name,standard)
        selection = group2selection(run,group)
        ns = length(selection)
        crunchers = Vector{NamedTuple}(undef,ns)
        for i in eachindex(selection)
            crunchers[i] = FCruncher(run[selection[i]],method,fit)
        end
        cruncher_groups[standard] = (anchor=anchor,crunchers=crunchers)
    end

    # initialise the parameters
    init = fill(0.0,method.ndrift)
    if (method.ndown>0)
        init = vcat(init,fill(0.0,method.ndown))
    end
    if isfinite(method.PAcutoff)
        init = vcat(init,fill(0.0,method.ndrift))
    end

    # define the objective function
    objective = (par) -> SS(par,method,cruncher_groups;
                            verbose=verbose)

    # fit the model
    optimum = Optim.optimize(objective,init)

    if verbose
        println("Drift and downhole fractionation correction:\n")
        println(optimum)
    else
        if optimum.stopped_by.time_limit
            @warn "Reached the maximum number of iterations " *
                "before achieving convergence. " *
                "Reduce the order of the polynomials or fix " *
                "the mass fractionation and try again."
        end
        if hasproperty(optimum.stopped_by,:ls_failed) &&
            optimum.stopped_by.ls_failed
            @warn "Least squares algorithm did not converge."
        end
    end

    # update the fit
    solution = Optim.minimizer(optimum)
    par2Gfit!(fit,method,solution)
    fractionation_error!(fit,objective,solution)
end
function fractionation_error!(fit::Gfit,
                              objective::Function,
                              solution::AbstractVector)
    H = FiniteDiff.finite_difference_hessian(objective, solution)
    if rank(H)==size(H,1)
        fit.covmat = inv(H/2)
    else
        fit.covmat = pinv(H/2)
    end
end

"""
     Cs  sum(S_i X_i)
[f = -- --------------]
     C    sum(S_i^2)

with C, Cs = reference concentrations of the elements and internal standard
     X_i, S_i = blank-corrected measurement for elements and internal standard
"""
function fractionation!(fit::Cfit,
                        method::Cmethod,
                        run::Vector{Sample};
                        kwargs...)
    num = fit.blank[1:1,:] .* 0.0
    den = copy(num)
    internal = method.internal[1]
    for standard in method.standards
        selection = getIndicesInGroup(run,standard)
        dats = [swinData(samp) for samp in run[selection]]
        for dat in dats
            bt = polyVal(fit.blank,dat.t)
            X = getSignals(dat) .- bt
            C = getConcentrations(method,standard)
            num[1,:] = Vector(num[1,:]) + sum.(eachcol(C[1,internal].*X.*X[:,internal]))
            den[1,:] = Vector(den[1,:]) + sum.(eachcol(C.*(X[:,internal].^2)))
        end
    end
    fit.par = num./den
end
export fractionation!

function par2Gfit!(fit::Gfit,
                   method::Gmethod,
                   par::AbstractVector)
    fit.drift = par[1:method.ndrift]
    fit.down = vcat(0.0,par[method.ndrift+1:method.ndrift+method.ndown])
    fit.adrift = isnothing(method.PAcutoff) ? fit.drift : par[end-method.ndrift+1:end]
end
function par2fit(par::AbstractVector,
                 method::Gmethod)
    fit = Gfit(method)
    par2Gfit!(fit,method,par)
    return fit
end

function FCruncher(samp::Sample,
                   method::Gmethod,
                   fit::Gfit)

    dat = swinData(samp)
    
    pm = dat[:,method.P.channel]
    Dm = dat[:,method.D.channel]
    bm = dat[:,method.d.channel]

    t = dat.t
    T = dat.T

    bpt = polyVal(fit.blank[:,method.P.channel],t)
    bDt = polyVal(fit.blank[:,method.D.channel],t)
    bbt = polyVal(fit.blank[:,method.d.channel],t)

    pmb = pm - bpt
    Dmb = Dm - bDt
    bmb = bm - bbt

    Ip = interference_correction(dat,method.P.interferences,fit.bias)
    ID = interference_correction(dat,method.D.interferences,fit.bias)
    Ib = interference_correction(dat,method.d.interferences,fit.bias)

    mf = bias_correction(t,fit.bias;num=method.d.proxy,den=method.D.proxy)

    sig = hcat(pmb,Dmb,bmb)
    covmat = df2cov(sig)
    vp = covmat[1,1]
    vD = covmat[2,2]
    vb = covmat[3,3]
    spD = covmat[1,2]
    spb = covmat[1,3]
    sDb = covmat[2,3]
    
    bd = iratio(method.d.proxy,method.d.ion)

    return (pmb=pmb,Dmb=Dmb,bmb=bmb,
            bpt=bpt,bDt=bDt,bbt=bbt,
            vp=vp,vD=vD,vb=vb,
            spD=spD,spb=spb,sDb=sDb,
            Ip=Ip,ID=ID,Ib=Ib,
            mf=mf,bd=bd,t=t,T=T)
    
end
export FCruncher