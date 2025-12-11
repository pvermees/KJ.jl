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
              covmat::AbstractMatrix = zeros(length([drift,down]),
                                             length([drift,down])))
    return Gfit(blank,drift,down,drift,covmat)
end

function Cfit()
    blank = DataFrame()
    par = DataFrame()
    return Cfit(blank,par)
end