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
    x = nothing
    y = nothing
    if add_xy && "x" in names(dat) && "y" in names(dat)
        x=dat.x
        y=dat.y
    end
    return (P=Phat, D=Dhat, d=dhat, x=x, y=y)
end
export atomic