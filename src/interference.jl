function interference_correction(dat::AbstractDataFrame,
                                 interferences::Set{AbstractInterference};
                                 bias::AbstractDict=Dict{String,Bias}(),
                                 blank::AbstractDataFrame=DataFrame())
    out = fill(0.0,size(dat,1))
    for interference in interferences
        out .+= interference_correction(dat,interference;
                                        bias=bias,blank=blank)
    end
    return out
end
function interference_correction(dat::AbstractDataFrame,
                                 interference::Interference;
                                 bias::AbstractDict=Dict{String,Bias}(),
                                 blank::AbstractDataFrame=DataFrame())
    ch = interference.channel
    meas = dat[:,ch]
    blk = polyVal(blank[:,ch],dat.t)
    y = iratio(interference.ion,interference.proxy)
    mf = bias4interference(dat,interference,bias)
    return @. (meas - blk) * y * mf
end
function interference_correction(dat::AbstractDataFrame,
                                 interference::REEInterference;
                                 bias::AbstractDict=Dict{String,REEBias}(),
                                 blank::AbstractDataFrame=DataFrame())
    ch = interference.proxy
    meas = dat[:,ch]
    blk = polyVal(blank[:,ch],dat.t)
    mf = polyFac(bias[ch].par,dat.t)
    return @. (meas - blk) * mf
end
export interference_correction

function bias4interference(dat::AbstractDataFrame,
                           interference::Interference,
                           bias::Dict{String,AbstractBias})
    element = channel2element(interference.ion)
    if haskey(bias,element)
        mf = bias_correction(bias[element],interference.ion,
                             interference.proxy,dat.t)
    else
        mf = fill(1.0,size(dat,1))
    end
    return mf
end
function bias4interference(dat::AbstractDataFrame,
                           interference::REEInterference,
                           bias::Dict{String,AbstractBias})
    key = interference.bias_key
    if key in names(bias)
        return polyFac(bias[:,key],dat.t)
    else
        error("No bias correction supplied to REEInterference.")
    end
end