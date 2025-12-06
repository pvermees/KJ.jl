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
                       mu_p::Real=0.5,
                       sigma_p::Real=0.2,
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
    FT = polyFac(fit.down,T)
    bt = polyVal(fit.blank,t)
    for col in eachcol(bt)
        col .+= (sqrt(Distributions.median(col)) .* randn(nsig+nblk))
    end
    p = fill(mu_p + randn() * sigma_p, nsig)
    x = p.*x0
    y = @. y0 + (y1-y0)*x/x0
    isig = nblk.+(1:nsig)
    P = x.*D
    d = y.*D
    ch = getChannels(method;as_tuple=true)
    Pm = bt[:,ch.P]
    Dm = bt[:,ch.D]
    dm = bt[:,ch.d]
    Dm[isig] .+= D
    Pm[isig] .+= P.*ft[isig].*FT[isig]
    dm[isig] .+= d
    eD = Distributions.median(Dm[isig]).*relerr_D.*randn(nsig)
    eP = Distributions.median(Pm[isig]).*relerr_P.*randn(nsig)
    ed = Distributions.median(dm[isig]).*relerr_d.*randn(nsig)
    Pm[isig] .+= eP
    Dm[isig] .+= eD
    dm[isig] .+= ed
    channels = getChannels(method;as_tuple=true)
    dat = DataFrame("Time [sec]" => spot_time,
                    channels.P => Pm,
                    channels.D => Dm,
                    channels.d => dm,
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
    a = method.anchors
    x0_std = a[first(keys(a))].x0
    y0_std = a[first(keys(a))].y0
    ch = getChannels(method;as_tuple=true)
    b = D/1000.0
    blank = DataFrame(ch.P => [b], ch.D => [b], ch.d => [b])
    covmat = zeros(length(drift)+length(down)-1,
                   length(drift)+length(down)-1)
    fit = Gfit(blank,drift,down,drift,covmat)
    run = [
        random_sample(method,fit;i=1,n=4,
                      nblk=nblk,nsig=nsig,sname="BP_1",
                      dtime=DateTime("2025-01-01T08:00:00"),
                      spot_time=spot_time,t0=t0,bwin=bwin,swin=swin,
                      mu_p=0.5,sigma_p=0.1,x0=x0_std,y0=y0_std,
                      D=D,relerr_D=relerr_D,relerr_P=relerr_P,relerr_d=relerr_d,
                      group="BP_gt"),
        random_sample(method,fit;i=2,n=4,
                      nblk=nblk,nsig=nsig,sname="Hog_1",
                      dtime=DateTime("2025-01-01T08:01:00"),
                      spot_time=spot_time,t0=t0,bwin=bwin,swin=swin,
                      mu_p=0.5,sigma_p=0.1,x0=x0_smp,y0=y0_smp,
                      D=D,relerr_D=relerr_D,relerr_P=relerr_P,relerr_d=relerr_d),
        random_sample(method,fit;i=3,n=4,
                      nblk=nblk,nsig=nsig,sname="NIST612_1",
                      dtime=DateTime("2025-01-01T08:02:00"),
                      spot_time=spot_time,t0=t0,bwin=bwin,swin=swin,
                      mu_p=0.5,sigma_p=0.0,y0=2*y0_glass,
                      D=D,relerr_D=relerr_D,relerr_P=relerr_P,relerr_d=relerr_d,
                      group="NIST612"),
        random_sample(method,fit;i=4,n=4,
                      nblk=nblk,nsig=nsig,sname="BP_2",
                      dtime=DateTime("2025-01-01T08:03:00"),
                      spot_time=spot_time,t0=t0,bwin=bwin,swin=swin,
                      mu_p=0.5,sigma_p=0.1,x0=x0_std,y0=y0_std,
                      D=2*D,relerr_D=relerr_D,relerr_P=relerr_P,relerr_d=relerr_d,
                      group="BP_gt")
    ]
    return run, fit
end