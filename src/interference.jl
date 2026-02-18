"""
    interference_correction(dat::AbstractDataFrame, interferences::AbstractDict; bias=Dict(), blank=DataFrame())
    interference_correction(dat::AbstractDataFrame, ion::AbstractString, interference::Interference; bias=Dict(), blank=DataFrame())
    interference_correction(dat::AbstractDataFrame, proxy_channel::AbstractString, interference::REEInterference; bias=Dict(), blank=DataFrame())

Apply interference corrections to measured data.

# Arguments
- `dat`: DataFrame containing measured data with a time column
- `interferences`: Dictionary of interferences to correct for, or a single ion/channel and interference
- `bias`: Optional dictionary of bias corrections
- `blank`: Optional DataFrame with blank corrections

# Returns
- Vector of interference-corrected values
"""
function interference_correction(dat::AbstractDataFrame,
                                 interferences::AbstractDict;
                                 bias::AbstractDict=Dict(),
                                 blank::AbstractDataFrame=DataFrame())
    out = fill(0.0,size(dat,1))
    for (key,interference) in interferences
        out .+= interference_correction(dat,key,interference;
                                        bias=bias,blank=blank)
    end
    return out
end
function interference_correction(dat::AbstractDataFrame,
                                 ion::AbstractString,
                                 interference::Interference;
                                 bias::AbstractDict=Dict(),
                                 blank::AbstractDataFrame=DataFrame())
    ch = interference.channel
    meas = dat[:,ch]
    blk = polyVal(blank[:,ch],dat.t)
    y = iratio(ion,interference.proxy)
    mf = bias4interference(dat,ion,interference,bias)
    return @. (meas - blk) * y * mf
end
function interference_correction(dat::AbstractDataFrame,
                                 proxy_channel::AbstractString,
                                 interference::REEInterference;
                                 bias::AbstractDict=Dict(),
                                 blank::AbstractDataFrame=DataFrame())
    meas = dat[:,proxy_channel]
    blk = polyVal(blank[:,proxy_channel],dat.t)
    mf = polyFac(bias[proxy_channel].par,dat.t)
    return @. (meas - blk) * mf
end
export interference_correction

function bias4interference(dat::AbstractDataFrame,
                           ion::AbstractString,
                           interference::Interference,
                           bias::AbstractDict)
    element = channel2element(ion)
    if haskey(bias,element)
        mf = bias_correction(bias[element],ion,
                             interference.proxy,dat.t)
    else
        mf = fill(1.0,size(dat,1))
    end
    return mf
end