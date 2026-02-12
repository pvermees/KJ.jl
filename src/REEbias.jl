function fit_bias(run::Vector{Sample},
                  method::Gmethod,
                  interference::REEInterference,
                  blank::AbstractDataFrame)
    crunchers = NamedTuple[]
    for group in interference.standards
        selection = group2selection(run,group)
        for i in eachindex(selection)
            samp = run[selection[i]]
            cruncher = BCruncher(samp,interference,blank)
            push!(crunchers,cruncher)
        end
    end

    init = [init_bias(crunchers);fill(0.0,method.nbias-1)]
    objective = (par) -> SS(par,crunchers)
    optimum = Optim.optimize(objective,init)
    fit = Optim.minimizer(optimum)
    return REEBias(fit)
end

function init_bias(crunchers::AbstractVector)
    Dsum = 0.0
    bsum = 0.0
    for cruncher in crunchers
        Dsum += sum(cruncher[:Dmb])
        bsum += sum(cruncher[:bmb])
    end
    return log(bsum) - log(Dsum)
end

function add_bias!(bias::AbstractDict,
                   run::Vector{Sample},
                   method::Gmethod,
                   blank::AbstractDataFrame,
                   bias_key::AbstractString,
                   interference::REEInterference)
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