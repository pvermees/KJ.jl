function KJfit(method::Gmethod)
    return Gfit(method)
end

function KJfit(method::Cmethod)
    return Cfit()
end

function init_bias(method::Gmethod)
    colnames = [collect(keys(method.fractionation.bias));
                collect(keys(method.interference.bias))]
    data = zeros(method.nbias,length(colnames))
    return DataFrame(data,colnames)
end

function Gfit(method::Gmethod;
              blank::AbstractDataFrame = DataFrame(fill(0.0,method.nblank,3),
                                                   getChannels(method)),
              drift::AbstractVector = fill(0.0,method.ndrift),
              down::AbstractVector = fill(0.0,method.ndown),
              adrift::AbstractVector = drift,
              covmat::AbstractMatrix = zeros(length([drift,down]),
                                             length([drift,down])),
              bias::AbstractDataFrame=init_bias(method))
    return Gfit(blank,drift,down,adrift,covmat,bias)
end

function Cfit()
    blank = DataFrame()
    par = DataFrame()
    return Cfit(blank,par)
end