function fit_bias(run::Vector{Sample},
                  method::Gmethod,
                  element::AbstractString,
                  blank::AbstractDataFrame)
    calibration = method.bias[element]
    m1 = get_proxy_isotope(calibration.num.ion,element)
    m2 = get_proxy_isotope(calibration.den.ion,element)
    cruncher_groups = Dict()
    for group in calibration.standards
        standard = method.groups[group]
        selection = group2selection(run,group)
        ns = length(selection)
        crunchers = Vector{NamedTuple}(undef,ns)
        for i in eachindex(selection)
            crunchers[i] = BCruncher(run[i],method,element,blank)
        end
        cruncher_groups[standard] = crunchers
    end
end
export fit_bias

function bias_correction(t::AbstractVector,
                         par::AbstractVector,
                         mass1::Int,
                         mass2::Int)
    beta = (1/mass1 - 1/(mass1+1))/(1/mass1 - 1/mass2)
    return polyFac(par,t).^beta
end

function BCruncher(samp::Sample,
                   method::Gmethod,
                   element::AbstractString,
                   blank::AbstractDataFrame)
    dat = swinData(samp)
    t = dat.T
    calibration = method.bias[element]
    Nm = dat[:,calibration.num.channel]
    Dm = dat[:,calibration.den.channel]
    return (t=t,Nm=Nm,Dm=Dm)
end