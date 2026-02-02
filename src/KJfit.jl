function KJfit(method::Gmethod)
    return Gfit(method)
end

function KJfit(method::Cmethod)
    return Cfit()
end

function Gfit(method::Gmethod;
              blank::AbstractDataFrame = DataFrame(fill(0.0,method.nblank,3),
                                                   getChannels(method)),
              drift::AbstractVector = fill(0.0,method.ndrift),
              down::AbstractVector = fill(0.0,method.ndown),
              adrift::AbstractVector = drift,
              covmat::AbstractMatrix = zeros(length([drift,down]),
                                             length([drift,down])),
              bias::AbstractDict = Dict{String,AbstractBias}())
    return Gfit(blank,drift,down,adrift,covmat,bias)
end

function Cfit()
    blank = DataFrame()
    par = DataFrame()
    return Cfit(blank,par)
end