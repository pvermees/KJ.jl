function KJfit(method::KJmethod)
    blank = DataFrame()
    drift = fill(0.0,method.ndrift)
    down = fill(0.0,method.ndown)
    return KJfit(blank,drift,down,drift)
end
