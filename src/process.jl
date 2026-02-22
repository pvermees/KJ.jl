"""
    process!(run::Vector{Sample}, method::KJmethod; reject_outliers=true, verbose=false)

Process a complete run of LA-ICP-MS data.

This function performs the full data reduction workflow:
1. Assign group labels to samples based on the method
2. Detect and flag outliers (if reject_outliers=true)
3. Initialize fit object
4. Fit blank corrections
5. Fit bias corrections (for Gmethod)
6. Fit drift and downhole fractionation corrections

# Arguments
- `run`: Vector of samples to process
- `method`: Geochronology or concentration method definition
- `reject_outliers`: Whether to automatically detect and flag outliers (default: true)
- `verbose`: Print detailed diagnostic information (default: false)

# Returns
- Fitted KJfit object containing all correction parameters
"""
function process!(run::Vector{Sample},
                  method::KJmethod;
                  setGroup::Bool=true,
                  reject_outliers::Bool=true,
                  verbose::Bool=false)
    if setGroup
        setGroup!(run,method)
    end
    if reject_outliers
        detect_outliers!(run,method)
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