function atomic(samp::Sample,
                method::Gmethod,
                fit::Gfit;
                add_xy::Bool=false)
    dat = swinData(samp;add_xy=add_xy)
    c = Cruncher(samp,method,fit)
    ft, FT = ft_FT(c,method,fit)
    Phat = @. c.pmb/(ft*FT)
    Dhat = @. c.Dombi
    dhat = @. c.bomb/c.bd
    x = nothing
    y = nothing
    if add_xy && "x" in names(dat) && "y" in names(dat)
        x=dat.x
        y=dat.y
    end
    return (P=Phat, D=Dhat, d=dhat, x=x, y=y)
end
export atomic
