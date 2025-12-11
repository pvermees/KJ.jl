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

function fitBlanks!(fit::KJfit,
                    method::KJmethod,
                    run::Vector{Sample})
    fit.blank = fitBlanks(run;
                          nblank=method.nblank)
end
export fitBlanks!
