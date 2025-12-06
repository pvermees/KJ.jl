function KJfit(method::Gmethod)
    return Gfit(method)
end

function KJfit(method::Cmethod)
    return Cfit()
end

function Gfit(method::Gmethod)
    blank = DataFrame()
    drift = fill(0.0,method.ndrift)
    down = fill(0.0,method.ndown)
    covmat = zeros(method.ndrift+method.ndown,
                   method.ndrift+method.ndown)
    return Gfit(blank,drift,down,drift,covmat)
end

function Cfit()
    blank = DataFrame()
    par = DataFrame()
    return Cfit(blank,par)
end