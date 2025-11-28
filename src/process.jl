function process!(run::Vector{Sample},
                  method::KJmethod,
                  standards::AbstractDict;
                  reject_outliers::Bool=true,
                  verbose::Bool=false)
    if reject_outliers
        ch = collect(values(getChannels(method)))
        detect_outliers!(run;channels=ch)
    end
    fit = KJfit(method)
    setStandards!(run,method;standards=standards)
    fitBlanks!(fit,method,run)
    fractionation!(fit,method,run;verbose=verbose)
    return fit
end

function process!(run::Vector{Sample},
                  internal::Tuple,
                  glass::AbstractDict;
                  nblank::Integer=2,
                  reject_outliers::Bool=true)
    setGroup!(run,glass)
    if reject_outliers
        detect_outliers!(run)
    end
    blank = fitBlanks(run;nblank=nblank)
    fit = fractionation(run,blank,internal,glass)
    return blank, fit
end
export process!
