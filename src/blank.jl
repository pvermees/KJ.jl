"""
    blank!(fit::KJfit, method::KJmethod, run::Vector{Sample})

Fit blank corrections for a run.

Extracts data from blank windows and fits polynomial background models
for each channel.

# Arguments
- `fit`: Fit object to populate with blank parameters
- `method`: Method definition
- `run`: Vector of samples
"""
function blank!(fit::KJfit,
                method::KJmethod,
                run::Vector{Sample})
    fit.blank = fitBlanks(run;
                          nblank=method.nblank)
end
export blank!

"""
    fitBlanks(run::Vector{Sample}; nblank=2)

Fit polynomial blank models to blank window data.

# Arguments
- `run`: Vector of samples
- `nblank`: Order of polynomial (default: 2 for quadratic)

# Returns
- DataFrame of polynomial coefficients for each channel
"""
function fitBlanks(run::Vector{Sample};
                   nblank=2)
    blk = reduce(vcat, bwinData(samp) for samp in run)
    channels = getChannels(run)
    nc = length(channels)
    bpar = DataFrame(zeros(nblank,nc),channels)
    good = .!blk.outlier
    for channel in channels
        bpar[:,channel] = polyFit(blk.t[good],blk[good,channel],nblank)
    end
    return bpar
end
export fitBlanks

function init_blank(method::KJmethod)
    channels = getChannels(method)
    nc = length(channels)
    return DataFrame(fill(0.0,method.nblank,nc), channels)
end
