function atomic(samp::Sample,
                method::KJmethod,
                fit::KJfit;
                add_xy::Bool=false)
    dat = swinData(samp;add_xy=add_xy)
    sig = getSignals(dat)
    covmat = df2cov(sig)
    c = Cruncher(samp,method,fit)
    ft, FT = ft_FT(c,method,fit)
    Phat = @. c.pmb/(ft*FT)
    Dhat = @. c.Dombi/mf
    dhat = @. c.bomb
    if add_xy
        if all([cname in names(dat) for cname in ["x";"y"]])
            return (Phat, Dhat, dhat, dat.x, dat.y)
        else
            return (Phat, Dhat, dhat, nothing, nothing)
        end
    else
        return (Phat, Dhat, dhat)
    end
end
export atomic
