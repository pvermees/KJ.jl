"""
fitBlanks(run::Vector{Sample};
          nblank=2)

Fit a polynomial to the blanks in run.
"""
function fitBlanks(run::Vector{Sample};
                   nblank=2)
    blks = pool(run;blank=true)
    blk = reduce(vcat,blks)
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
