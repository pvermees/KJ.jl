function fractionation!(fit::Gfit,
                        method::Gmethod,
                        run::Vector{Sample};
                        verbose::Bool=false)

    # extract the grouped data for the SS function from the run
    cruncher_groups = Dict()
    for refmat in keys(method.anchors)
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
function fractionation!(fit::Cfit,
                        method::Cmethod,
                        run::Vector{Sample};
                        kwargs...)
    ne = ncol(method.elements)
    num = den = fill(0.0,ne-1)
    internal = method.internal[1]
    for SRM in method.refmats
        selection = getIndicesInGroup(run,SRM)
        dats = [swinData(samp) for samp in run[selection]]
        concs = method.concentrations
        for dat in dats
            bt = polyVal(fit.blank,dat.t)
            sig = getSignals(dat)
            Xm = sig[:,Not(internal)]
            Sm = sig[:,internal]
            bXt = bt[:,Not(internal)]
            bSt = bt[:,internal]
            S = Sm.-bSt
            Cint = concs[:,internal]
            Coth = concs[:,Not(internal)]
            R = collect((Coth./Cint)[1,:])
            num += sum.(eachcol(R'.*(Xm.-bXt).*S))
            den += sum.(eachcol((R'.*S).^2))
        end
    end
    fit.par = fit.blank[1:1,Not(internal)]
    fit.par[1,:] = num./den
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