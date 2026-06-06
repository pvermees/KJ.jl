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
                   nblank::Int = 2)
    if nblank < 1
        return fitPiecewiseBlanks(run)
    else
        return fitPolyBlanks(run, nblank)
    end
end
export fitBlanks

function fitPiecewiseBlanks(run::Vector{Sample})
    ns = length(run)
    channels = getChannels(run)
    out = hcat(DataFrame(:sample => fill("",ns)),
               DataFrame([name => zeros(ns) for name in channels]))
    for i in eachindex(run)
        samp = run[i]
        blk = bwinData(samp)
        out[i,:sample] = samp.sname
        out[i,channels] = Statistics.mean.(eachcol(blk[:,channels]))
    end
    return out
end

function fitPolyBlanks(run::Vector{Sample},nblank::Int)
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

function isPolyBlank(blank::AbstractDataFrame)
    return names(blank)[1] !== "sample"
end

function init_blank(method::KJmethod)
    channels = getChannels(method)
    nc = length(channels)
    return DataFrame(zeros(method.nblank,nc), channels)
end

function plot(blk::AbstractDataFrame,
              run::Vector{Sample};
              plot_options...)
    ns = length(run)
    t, y, yf, yl, yu = zeros(ns), zeros(ns), zeros(ns), zeros(ns), zeros(ns)
    polyBlank = isPolyBlank(blk)
    for i in eachindex(run)
        df = bwinData(run[i])
        sig = getSignals(df)
        if polyBlank
            fit = polyVal(blk,df.t)
            yf[i] = Statistics.mean(sum.(eachrow(fit)))
        else
            yf[i] = sum(blk[i,Not(:sample)])
        end
        t[i] = Statistics.median(df.t)
        y[i] = Statistics.mean(sum.(eachrow(sig)))
        yl[i] = Statistics.quantile(sum.(eachrow(sig)), 0.025)
        yu[i] = Statistics.quantile(sum.(eachrow(sig)), 0.975)
    end
    p = Plots.scatter(t, y, yerror=(y-yl, yu-y), marker=:circle, label=false; plot_options...)
    o = detect_outliers(y)
    annotate!(p, t[o], yu[o], text.(o,7,:center,:bottom))
    Plots.plot!(p, t, yf, label=nothing, color=:red)
    return p
end