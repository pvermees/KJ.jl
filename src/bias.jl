function bias(run::Vector{Sample},
              method::Gmethod,
              blank::AbstractDataFrame)
    for (target,interferences) in method.interference.ions
        for interference in interferences

        end
    end
    return nothing
end

function bias(run::Vector{Sample},
              fractionation::Fractionation)
    c = BCruncher(run,fractionation)
     return nothing
end
export bias

function BCruncher(run::Vector{Sample},
                   method::Gmethod,
                   blank::AbstractDataFrame)
    return nothing
end

function BCruncher(run::Vector{Sample},
                   fractionation::Fractionation)
    return nothing
end