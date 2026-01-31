function interference_correction(dat::AbstractDataFrame,
                                 interferences::Set{AbstractInterference};
                                 bias::AbstractDataFrame=DataFrame(),
                                 blank::AbstractDataFrame=DataFrame())
    out = fill(0.0,size(dat,1))
    for interference in interferences
        out .+= interference_correction(dat,interference;
                                        bias=bias,blank=blank)
    end
    return out
end
function interference_correction(dat::AbstractDataFrame,
                                 interference::AbstractInterference;
                                 bias::AbstractDataFrame=DataFrame(),
                                 blank::AbstractDataFrame=DataFrame())
    meas = dat[:,interference.channel]
    blk = polyVal(blank[:,interference.channel],dat.t)
    y = iratio(interference.ion,interference.proxy)
    mf = bias4interference(dat,interference,bias)
    return @. (meas - blk) * y * mf
end
export interference_correction

function bias4interference(dat::AbstractDataFrame,
                           interference::Interference,
                           bias::AbstractDataFrame)
    element = channel2element(interference.ion)
    if element in names(bias)
        mf = bias_correction(bias[:,element],interference.ion,
                             interference.proxy,dat.t)
    else
        mf = fill(1.0,size(dat,1))
    end
    return mf
end
function bias4interference(dat::AbstractDataFrame,
                           interference::REEInterference,
                           bias::AbstractDataFrame)
    key = interference.bias_key
    if key in names(bias)
        return polyFac(bias[:,key],dat.t)
    else
        error("No bias correction supplied to REEInterference.")
    end
end