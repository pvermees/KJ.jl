"""
    fractionation!(fit::Gfit, method::Gmethod, run::Vector{Sample}; verbose=false)
    fractionation!(fit::Cfit, method::Cmethod, run::Vector{Sample}; kwargs...)

Fit drift and downhole fractionation corrections.

For geochronology methods (Gmethod), fits polynomial time-dependent drift and
downhole corrections using reference materials. For concentration methods (Cmethod),
fits sensitivity factors using internal standards.

# Arguments
- `fit`: Fit object to populate with fractionation parameters
- `method`: Method definition
- `run`: Vector of samples
- `verbose`: Print detailed optimization information (default: false)
"""
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
    init = zeros(method.ndrift)
    if (method.ndown>0)
        init = vcat(init,zeros(method.ndown))
    end
    if isfinite(method.PAcutoff)
        init = vcat(init,zeros(method.ndrift))
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
    channels = getChannels(run)
    num = DataFrame(zeros(1, length(channels)), channels)
    den = DataFrame(zeros(1, length(channels)), channels)
    internal = method.internal[1]
    for (group,standard) in method.groups
        selection = getIndicesInGroup(run,group)
        for samp in run[selection]
            dat = swinData(samp)
            bt = predict(samp,fit.blank;t=dat.t)
            X = getSignals(dat) .- bt
            S = X[:,internal]
            C = getConcentrations(method,standard)
            Cs = C[1,internal]
            num[1,:] = Vector(num[1,:]) + sum.(eachcol(Cs.*X.*S))
            den[1,:] = Vector(den[1,:]) + sum.(eachcol(C.*(S.^2)))
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
    fit.adrift = isfinite(method.PAcutoff) ? par[end-method.ndrift+1:end] : fit.drift
end
function par2fit(par::AbstractVector,
                 method::Gmethod)
    fit = Gfit(method)
    par2Gfit!(fit,method,par)
    return fit
end

"""
    FCruncher(samp::Sample, method::Gmethod, fit::Gfit)

Prepare per-sample arrays and covariance terms used by fractionation and
prediction routines.

# Returns
- Named tuple containing blank-corrected intensities, covariance terms,
  interference corrections, and time vectors.
"""
function FCruncher(samp::Sample,
                   method::Gmethod,
                   fit::Gfit)

    dat = swinData(samp)
    
    pm = dat[:,method.P.channel]
    Dm = dat[:,method.D.channel]
    bm = dat[:,method.d.channel]

    t = dat.t
    T = dat.T

    blk = predict(samp,fit.blank;t=t)
    bpt = blk[:,method.P.channel]
    bDt = blk[:,method.D.channel]
    bbt = blk[:,method.d.channel]

    pmb = pm - bpt
    Dmb = Dm - bDt
    bmb = bm - bbt

    Ip = interference_correction(dat,method.P.interferences;
                                 bias=fit.bias,blank=fit.blank)
    ID = interference_correction(dat,method.D.interferences;
                                 bias=fit.bias,blank=fit.blank)
    Ib = interference_correction(dat,method.d.interferences,
                                 bias=fit.bias,blank=fit.blank)

    Delement = channel2element(method.D.ion)
    if haskey(fit.bias,Delement)
        mf = bias_correction(fit.bias[Delement],
                             method.d.ion,method.D.ion,t)
    else
        mf = ones(length(t))
    end
    sig = hcat(pmb,Dmb,bmb)
    covmat = df2cov(sig)
    vp = covmat[1,1]
    vD = covmat[2,2]
    vb = covmat[3,3]
    spD = covmat[1,2]
    spb = covmat[1,3]
    sDb = covmat[2,3]
    
    bd = iratio(method.d.proxy,method.d.ion)
    if isnothing(bd)
        bd = 1.0
    end

    return (pmb=pmb,Dmb=Dmb,bmb=bmb,
            bpt=bpt,bDt=bDt,bbt=bbt,
            vp=vp,vD=vD,vb=vb,
            spD=spD,spb=spb,sDb=sDb,
            Ip=Ip,ID=ID,Ib=Ib,
            mf=mf,bd=bd,t=t,T=T)
    
end
export FCruncher