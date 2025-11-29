function KJfit(method::Gmethod)
    blank = DataFrame()
    drift = fill(0.0,method.ndrift)
    down = fill(0.0,method.ndown)
    return Gfit(blank,drift,down,drift)
end

function KJfit(method::Cmethod)
    blank = DataFrame()
    par = DataFrame()
    return Cfit(blank,par)
end