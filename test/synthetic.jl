function random_sample(method::Gmethod,
                       fit::Gfit;
                       i::Int,
                       n::Int,
                       nblk::Int,
                       nsig::Int,
                       sname::AbstractString,
                       dtime::DateTime,
                       spot_time::AbstractVector,
                       t0::Real,
                       bwin::AbstractArray,
                       swin::AbstractArray,
                       mu::Real=0.5,
                       sigma::Real=0.2,
                       x0::Real=1.0,
                       y0::Real=1.0,
                       y1::Real=0.0,
                       D::Real=1e4,
                       relerr_D::Real=0.1,
                       relerr_P::Real=0.1,
                       relerr_d::Real=0.1,
                       group::String="sample")
    Random.seed!(1)
    t = range((i-1)/n,i/n,length=nblk+nsig)
    T = (spot_time .- t0)./60
    ft = polyFac(fit.drift,t)
    hT = polyFac(fit.down,T)
    bt = polyVal(fit.blank,t)
    for col in eachcol(bt)
        col .+= (sqrt(Distributions.median(col)) .* randn(nsig+nblk))
    end
    px = fill(mu + randn() * sigma, nsig)
    x = px.*x0
    y = @. y0 + (y1-y0)*x/x0
    isig = nblk.+(1:nsig)
    P = x.*D
    d = y.*D
    pm = bt[:,method.P.channel]
    Dm = bt[:,method.D.channel]
    bm = bt[:,method.d.channel]
    pm[isig] .+= P.*ft[isig].*hT[isig]
    Dm[isig] .+= D
    bm[isig] .+= d*iratio(method.d.proxy,method.d.ion)
    ep = Distributions.median(pm[isig]).*relerr_P.*randn(nsig)
    eD = Distributions.median(Dm[isig]).*relerr_D.*randn(nsig)
    eb = Distributions.median(bm[isig]).*relerr_d.*randn(nsig)
    pm[isig] .+= ep
    Dm[isig] .+= eD
    bm[isig] .+= eb
    dat = DataFrame("Time [sec]" => spot_time,
                    method.P.channel => pm,
                    method.D.channel => Dm,
                    method.d.channel => bm,
                    "outlier" => falses(nblk+nsig),
                    "t" => t)
    return Sample(sname,dtime,dat,t0,bwin,swin,group)
end

function synthetic!(method::Gmethod;
                    drift::AbstractVector=[0.0],
                    down::AbstractVector=[0.0,0.0],
                    lambda::T,
                    t_std::T,
                    t_smp::T,
                    y0_smp::T,
                    y0_glass::T,
                    D::T=1e4,
                    relerr_D::T=0.1,
                    relerr_P::T=0.1,
                    relerr_d::T=0.1) where T <: Real
    nblk = 10
    nsig = 50
    spot_time = range(start=0.0,stop=60.0,length=nblk+nsig)
    t0 = 10.0
    bwin = [(1,9)]
    swin = [(11,59)]
    x0_std = 1/(exp(lambda*t_std)-1)
    x0_smp = 1/(exp(lambda*t_smp)-1)
    a = getAnchor(method.name,collect(method.standards)[1])
    y0_std = a.y0
    fit = Gfit(method;drift=drift,down=down)
    fit.blank[:,:] .= D./1000
    run = [
        random_sample(method,fit;i=1,n=4,
                      nblk=nblk,nsig=nsig,sname="BP_1",
                      dtime=DateTime("2025-01-01T08:00:00"),
                      spot_time=spot_time,t0=t0,bwin=bwin,swin=swin,
                      mu=0.5,sigma=0.1,x0=x0_std,y0=y0_std,
                      D=D,relerr_D=relerr_D,relerr_P=relerr_P,relerr_d=relerr_d,
                      group="BP"),
        random_sample(method,fit;i=2,n=4,
                      nblk=nblk,nsig=nsig,sname="Hog_1",
                      dtime=DateTime("2025-01-01T08:01:00"),
                      spot_time=spot_time,t0=t0,bwin=bwin,swin=swin,
                      mu=0.5,sigma=0.1,x0=x0_smp,y0=y0_smp,
                      D=D,relerr_D=relerr_D,relerr_P=relerr_P,relerr_d=relerr_d),
        random_sample(method,fit;i=3,n=4,
                      nblk=nblk,nsig=nsig,sname="NIST612_1",
                      dtime=DateTime("2025-01-01T08:02:00"),
                      spot_time=spot_time,t0=t0,bwin=bwin,swin=swin,
                      mu=0.5,sigma=0.0,y0=2*y0_glass,
                      D=D,relerr_D=relerr_D,relerr_P=relerr_P,relerr_d=relerr_d,
                      group="NIST612"),
        random_sample(method,fit;i=4,n=4,
                      nblk=nblk,nsig=nsig,sname="BP_2",
                      dtime=DateTime("2025-01-01T08:03:00"),
                      spot_time=spot_time,t0=t0,bwin=bwin,swin=swin,
                      mu=0.5,sigma=0.1,x0=x0_std,y0=y0_std,
                      D=D*2,relerr_D=relerr_D,relerr_P=relerr_P,relerr_d=relerr_d,
                      group="BP")
    ]
    return run, fit
end