function fractionation!(fit::Gfit,
                        method::Gmethod,
                        run::Vector{Sample};
                        verbose::Bool=false)

    # extract the grouped data for the SS function from the run
    cruncher_groups = Dict()
    for refmat in keys(method.standards)
        anchor = getAnchor(method.name,refmat)
        selection = group2selection(run,refmat)
        ns = length(selection)
        crunchers = Vector{Cruncher}(undef,ns)
        for i in eachindex(selection)
            crunchers[i] = Cruncher(run[selection[i]],method,fit)
        end
        cruncher_groups[refmat] = (anchor=anchor,crunchers=crunchers)
    end

    # initialise the parameters
    init = fill(0.0,method.ndrift)
    if (method.ndown>0)
        init = vcat(init,fill(0.0,method.ndown))
    end
    if !isnothing(method.PAcutoff)
        init = vcat(init,fill(0.0,method.ndrift))
    end

    # define the objective function
    objective = (par) -> SS(par,method,cruncher_groups;
                            verbose=verbose)

    # fit the model
    out = Optim.optimize(objective,init)
    if verbose
        println("Drift and downhole fractionation correction:\n")
        println(out)
    else
        if out.stopped_by.time_limit
            @warn "Reached the maximum number of iterations " *
                "before achieving convergence. " *
                "Reduce the order of the polynomials or fix " *
                "the mass fractionation and try again."
        end
        if hasproperty(out.stopped_by,:ls_failed) &&
            out.stopped_by.ls_failed
            @warn "Least squares algorithm did not converge."
        end
    end

    # update the fit
    par = Optim.minimizer(out)
    par2Gfit!(fit,par,method)

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
    for SRM in method.refmats
        selection = getIndicesInGroup(run,SRM)
        dats = [swinData(samp) for samp in run[selection]]
        for dat in dats
            bt = polyVal(fit.blank,dat.t)
            X = getSignals(dat) .- bt
            C = getConcentrations(method,SRM)
            num[1,:] = Vector(num[1,:]) + sum.(eachcol(C[1,internal].*X.*X[:,internal]))
            den[1,:] = Vector(den[1,:]) + sum.(eachcol(C.*(X[:,internal].^2)))
        end
    end
    fit.par = num./den
end
export fractionation!

function par2Gfit!(fit::Gfit,
                  par::AbstractVector,
                  method::KJmethod)
    fit.drift = par[1:method.ndrift]
    fit.down = vcat(0.0,par[method.ndrift+1:method.ndrift+method.ndown])
    fit.adrift = isnothing(method.PAcutoff) ? fit.drift : par[end-method.ndrift+1:end]
end
function par2fit(par::AbstractVector,
                 method::KJmethod)
    fit = KJfit(method)
    par2Gfit!(fit,par,method)
    return fit
end