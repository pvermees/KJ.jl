function process!(run::Vector{Sample},
                  method::KJmethod;
                  reject_outliers::Bool=true,
                  verbose::Bool=false)
    setGroup!(run,collect(keys(method.groups)))
    if reject_outliers
        ch = getChannels(method)
        detect_outliers!(run;channels=ch)
    end
    fit = KJfit(method)
    blank!(fit,method,run)
    if method isa Gmethod
        fit.bias = fit_bias(run,method,fit.blank)
    end
    fractionation!(fit,method,run;verbose=verbose)
    return fit
end
export process!