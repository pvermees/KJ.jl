function fit_bias(run::Vector{Sample},
                  method::Gmethod,
                  interference::REEInterference,
                  blank::AbstractDataFrame)
    cruncher_groups = Dict()
    for group in interference.standards
        standard = method.groups[group]
        selection = group2selection(run,group)
        ns = length(selection)
        crunchers = Vector{NamedTuple}(undef,ns)
        for i in eachindex(selection)
            samp = run[selection[i]]
            crunchers[i] = BCruncher(samp,interference,blank)
        end
        cruncher_groups[standard] = crunchers
    end

    init = [init_bias(cruncher_groups);fill(0.0,method.nbias-1)]
    objective = (par) -> SS(par,cruncher_groups)
    optimum = Optim.optimize(objective,init)
    fit = Optim.minimizer(optimum)
    return REEBias(fit)
end

function init_bias(cruncher_groups::AbstractDict)
    Dsum = 0.0
    bsum = 0.0
    for cruncher_group in values(cruncher_groups)
        for cruncher in cruncher_group
            Dsum += sum(cruncher[:Dmb])
            bsum += sum(cruncher[:bmb])
        end
    end
    return log(bsum) - log(Dsum)
end

function add_bias!(bias::AbstractDict,
                   run::Vector{Sample},
                   method::Gmethod,
                   blank::AbstractDataFrame,
                   interference::REEInterference)
    bias_key = interference.proxy
    bias[bias_key] = fit_bias(run,method,interference,blank)
end

function bias_correction(bias::REEBias;
                         t::AbstractVector,
                         other...)
    return polyFac(bias.par,t)
end

function BCruncher(samp::Sample,
                   interference::REEInterference,
                   blank::AbstractDataFrame)
    Dch = interference.REE
    bch = interference.REEO
    return BCruncher(samp,Dch,bch,blank)
end