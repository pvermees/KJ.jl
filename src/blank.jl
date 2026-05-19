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

function plot(blk::AbstractDataFrame,
              run::Vector{Sample})
    ns = length(run)
    t = fill(0.0,ns)
    y = fill(0.0,ns)
    yf = fill(0.0,ns)
    y025 = fill(0.0,ns)
    y975 = fill(0.0,ns)
    for i in eachindex(run)
        df = bwinData(run[i])
        sig = getSignals(df)
        fit = polyVal(blk,df.t)
        t[i] = Statistics.median(df.t)
        y[i] = Statistics.mean(sum.(eachrow(sig)))
        y025[i] = Statistics.quantile(sum.(eachrow(sig)), 0.025)
        y975[i] = Statistics.quantile(sum.(eachrow(sig)), 0.975)
        yf[i] = Statistics.mean(sum.(eachrow(fit)))
    end
    p = Plots.scatter(t,y, yerror=(y-y025, y975-y), marker=:circle, label=false)
    Plots.plot!(t,yf, label=false)
    p_top = twiny(p)
    plot!(
        p_top,
        xlims = xlims(p),
        xticks = (t,1:ns),
        label = false
    )
    return p
end