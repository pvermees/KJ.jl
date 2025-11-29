function process!(run::Vector{Sample},
                  method::KJmethod,
                  standards::AbstractDict;
                  reject_outliers::Bool=true,
                  verbose::Bool=false)
    setStandards!(run,method;standards=standards)
    if reject_outliers
        ch = collect(values(getChannels(method)))
        detect_outliers!(run;channels=ch)
    end
    fit = KJfit(method)
    fitBlanks!(fit,method,run)
    fractionation!(fit,method,run;verbose=verbose)
    return fit
end
export process!