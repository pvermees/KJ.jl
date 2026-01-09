function fractionation!(fit::Gfit,
                        method::Gmethod,
                        run::Vector{Sample};
                        verbose::Bool=false)

    # extract the grouped data for the SS function from the run
    cruncher_groups = Dict()
    for standard in method.fractionation.standards
        anchor = getAnchor(method.name,standard)
        selection = group2selection(run,standard)
        ns = length(selection)
        crunchers = Vector{FCruncher}(undef,ns)
        for i in eachindex(selection)
            crunchers[i] = FCruncher(run[selection[i]],
                                    method.fractionation,
                                    fit.blank)
        end
        cruncher_groups[standard] = (anchor=anchor,crunchers=crunchers)
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