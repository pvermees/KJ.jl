function interference_correction(dat::AbstractDataFrame,
                                 interferences::Set{AbstractInterference},
                                 bias::AbstractDataFrame)
    out = fill(0.0,size(dat,1))
    for interference in interferences
        out .+= interference_correction(dat,interference,bias)
    end
    return out
end
export interference_correction

function interference_correction(dat::AbstractDataFrame,
                                 interference::AbstractInterference,
                                 bias::AbstractDataFrame)
    # TODO
    return out = fill(0.0,size(dat,1))
end