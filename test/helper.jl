using Dates

function random_sample(i::Int,
                       n::Int;
                       nblk::Int,
                       nsig::Int,
                       sname::AbstractString,
                       dtime::DateTime,
                       spot_time::AbstractVector,
                       t0::Real,
                       bwin::AbstractArray,
                       swin::AbstractArray,
                       blank::AbstractDataFrame,
                       x0::Real=1.0,
                       y0::Real=1.0,
                       y1::Real=0.0,
                       mu_p::Real=0.5,
                       sigma_p::Real=0.2,
                       D::Real=1e4,
                       relerr_D::Real=0.1,
                       relerr_P::Real=0.1,
                       relerr_d::Real=0.1,
                       mfrac::Real=0.0,
                       drift::AbstractVector=[0.0],
                       down::AbstractVector=[0.0],
                       group::AbstractString="sample",
                       channels::AbstractDict)
    t = range((i-1)/n,i/n,length=nblk+nsig)
    T = (spot_time .- t0)./60
    mf = exp(mfrac)
    ft = polyFac(drift,t)
    FT = polyFac(down,T)
    bt = polyVal(blank,t)
    for col in eachcol(bt)
        col .+= (sqrt(median(col)) .* randn(nsig+nblk))
    end
    p = fill(mu_p + randn() * sigma_p, nsig)
    x = p.*x0
    y = @. y0 + (y1-y0)*x/x0
    isig = nblk.+(1:nsig)
    P = x.*D
    d = y.*D
    Dm = bt[:,"D"]
    Pm = bt[:,"P"]
    dm = bt[:,"d"]
    Dm[isig] .+= D.*mf
    Pm[isig] .+= P.*ft[isig].*FT[isig]
    dm[isig] .+= d
    eD = median(Dm[isig]).*relerr_D.*randn(nsig)
    eP = median(Pm[isig]).*relerr_P.*randn(nsig)
    ed = median(dm[isig]).*relerr_d.*randn(nsig)
    Pm[isig] .+= eP
    Dm[isig] .+= eD
    dm[isig] .+= ed
    dat = DataFrame("Time [sec]" => spot_time,
                    channels["P"] => Pm,
                    channels["D"] => Dm,
                    channels["d"] => dm,
                    "t" => t)
    return Sample(sname,dtime,dat,t0,bwin,swin,group)
end

function synthetic(;
                   lambda::T,
                   t_std::T,
                   y0_std::T,
                   t_smp::T,
                   y0_smp::T,
                   y0_glass::T,
                   channels::AbstractDict) where T <: Real
    nblk = 10
    nsig = 50
    spot_time = range(start=0.0,stop=60.0,length=nblk+nsig)
    t0 = 10.0
    bwin = [(1,9)]
    swin = [(12,59)]
    mfrac = 0.0
    drift = [0.0]
    down = [0.0]
    x0_std = 1/(exp(lambda*t_std)-1)
    x0_smp = 1/(exp(lambda*t_smp)-1)
    D_blank = 1.0
    blank = DataFrame("P" => [x0_std*D_blank*exp(down[1])*exp(drift[1]),0.1],
                      "D" => [D_blank*exp(mfrac[1]),0.1],
                      "d" => [y0_std*D_blank,0.1])
    run = [
        random_sample(1,4;
                      nblk=nblk,nsig=nsig,sname="BP_1",
                      dtime=DateTime("2025-01-01T08:00:00"),
                      spot_time=spot_time,t0=t0,bwin=bwin,swin=swin,blank=blank,
                      mu_p=0.5,sigma_p=0.1,x0=x0_std,y0=y0_std,
                      D=1e4,channels=channels,
                      mfrac=mfrac,drift=drift,down=down),
        random_sample(2,4;
                      nblk=nblk,nsig=nsig,sname="Hog_1",
                      dtime=DateTime("2025-01-01T08:01:00"),
                      spot_time=spot_time,t0=t0,bwin=bwin,swin=swin,blank=blank,
                      mu_p=0.5,sigma_p=0.1,x0=x0_smp,y0=y0_smp,
                      D=1e4,channels=channels,mfrac=mfrac,drift=drift,down=down),
        random_sample(3,4;
                      nblk=nblk,nsig=nsig,sname="NIST612_1",
                      dtime=DateTime("2025-01-01T08:02:00"),
                      spot_time=spot_time,t0=t0,bwin=bwin,swin=swin,blank=blank,
                      mu_p=0.5,sigma_p=0.1,y0=y0_glass,channels=channels,
                      D=1e4,relerr_D=0.1,relerr_P=0.01,relerr_d=0.1,
                      mfrac=mfrac,drift=drift,down=down),
        random_sample(4,4;
                      nblk=nblk,nsig=nsig,sname="BP_2",
                      dtime=DateTime("2025-01-01T08:03:00"),
                      spot_time=spot_time,t0=t0,bwin=bwin,swin=swin,blank=blank,
                      mu_p=0.5,sigma_p=0.1,x0=x0_std,y0=y0_std,
                      D=2e4,channels=channels,
                      mfrac=mfrac,drift=drift,down=down)
    ]
    return run, channels
end
