function interference_correction(dat::AbstractDataFrame,
                                 interferences::Set{AbstractInterference};
                                 blank::AbstractDataFrame=DataFrame(),
                                 bias::AbstractDataFrame=DataFrame())
    return fill(1.0,size(dat,1))
end
function interference_correction(dat::AbstractDataFrame,
                                 interference::Interference;
                                 blank::AbstractDataFrame=DataFrame(),
                                 bias::AbstractDataFrame=DataFrame())
    meas = dat[:,interference.channel]
    blk = polyVal(blank[:,interference.channel],dat.t)
    y = iratio(interference.ion,interference.proxy)
    element = channel2element(interference.ion)
    if element in names(bias)
        m1 = get_proxy_isotope(interference.ion)
        m2 = get_proxy_isotope(interference.proxy)
        mf = bias_correction(bias[:,element],m1,m2;dat.t)
    else
        mf = fill(1.0,length(meas))
    end
    return @. (meas - blk) * y * mf
end
export interference_correction