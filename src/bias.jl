function bias(run::Vector{Sample},
              interference::Interference;
              nbias::Int=1)
    return nothing
end

function bias(run::Vector{Sample},
              fractionation::Fractionation;
              nbias::Int=1)
    return nothing
end
export bias