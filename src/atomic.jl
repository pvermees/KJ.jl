"""
    atomic(samp::Sample, method::Gmethod, fit::Gfit; add_xy=false)

Calculate atomic abundances (P, D, d) from measured signals.

Applies all corrections (blank, interference, bias, fractionation) to convert
raw signal intensities to corrected atomic abundances.

# Arguments
- `samp`: Sample to process
- `method`: Geochronology method definition
- `fit`: Fitted correction parameters
- `add_xy`: Include x,y spatial coordinates if available (default: false)

# Returns
- DataFrame with fields P (parent), D (daughter), d (sister/normalizing isotope),
  and optionally x, y coordinates
"""
function atomic(samp::Sample,
                method::Gmethod,
                fit::Gfit;
                add_xy::Bool=false)
    dat = swinData(samp;add_xy=add_xy)
    c = FCruncher(samp,method,fit)
    ft, hT = ft_hT(fit,method.PAcutoff;c...)
    Phat = @. (c.pmb-c.Ip)/(ft*hT)
    Dhat = @. (c.Dmb-c.ID)
    dhat = @. (c.bmb-c.Ib)/(c.bd*c.mf)
    out = DataFrame(P=Phat, D=Dhat, d=dhat)
    if add_xy && "x" in names(dat) && "y" in names(dat)
        out.x = dat.x
        out.y = dat.y
    end
    return out
end
export atomic