using Dates

function random_sample(i::Int,
                       n::Int;
                       nblk::Int,
                       nsig::Int,
                       sname::AbstractString,
                       dtime::DateTime,
                       Time::AbstractVector,
                       t0::Real,
                       bwin::AbstractArray,
                       swin::AbstractArray,
                       blank::AbstractDataFrame,
                       x0::Real=1.0,
                       y0::Real=1.0,
                       y1::Real=0.0,
                       p::Real,
                       mu_D::Real,
                       sigma_D::Real,
                       mfrac::Real,
                       drift::AbstractVector=[0.0],
                       down::AbstractVector=[0.0],
                       group::AbstractString="sample")
    t = range((i-1)/n,i/n,length=nblk+nsig)
    T = (Time .- t0)./60
    mf = exp(mfrac)
    ismineral = true
    ft = polyFac(drift,t)
    FT = polyFac(down,T)
    bt = polyVal(blank,t)
    x = p*x0
    y = y0 + (y1-y0)*x/x0
    isig = nblk.+(1:nsig)
    D = mu_D .+ randn(nsig) .* sigma_D
    P = x.*D
    d = y.*D
    Dm = bt[:,"D"]
    Pm = bt[:,"P"]
    dm = bt[:,"d"]
    Dm[isig] .+= D.*mf
    Pm[isig] .+= P.*ft[isig].*FT[isig]
    dm[isig] .+= d
    dat = DataFrame("Time [sec]" => Time,
                    "Lu175 -> 175" => Pm,
                    "Hf176 -> 258" => Dm,
                    "Hf178 -> 260" => dm,
                    "t" => t)
    return Sample(sname,dtime,dat,t0,bwin,swin,group)
end

function synthetic(;lambda::T,t_std::T,y0_std::T,
                   t_smp::T,y0_smp::T,y0_glass::T) where T <: Real
    nblk = 10
    nsig = 50
    P = [0.0]
    D = [0.0]
    d = [0.0]
    Time = range(start=0.0,stop=60.0,length=nblk+nsig)
    blank = DataFrame("P" => [0.0,0.1],
                      "D" => [0.0,0.1],
                      "d" => [0.0,0.1])
    Time = range(0.0,60.0,length=nblk+nsig)
    t0 = 10.0
    bwin = [(1,9)]
    swin = [(11,59)]
    mfrac = 0.0
    drift = [0.0]
    down = [0.0]
    x0_std = log(lambda*t_std)-1
    x0_smp = log(lambda*t_smp)-1
    run = [
        random_sample(1,4;
                      nblk=nblk,nsig=nsig,sname="BP_1",
                      dtime=DateTime("2025-01-01T08:00:00"),
                      Time=Time,t0=t0,bwin=bwin,swin=swin,blank=blank,
                      mu_D=1000,sigma_D=50,x0=x0_std,y0=y0_std,
                      p=0.5,mfrac=mfrac,drift=drift,down=down,
                      group="BP_gt"),
        random_sample(2,4;
                      nblk=nblk,nsig=nsig,sname="Hog_1",
                      dtime=DateTime("2025-01-01T08:01:00"),
                      Time=Time,t0=t0,bwin=bwin,swin=swin,blank=blank,
                      mu_D=2000,sigma_D=75,x0=x0_smp,y0=y0_smp,
                      p=0.5,mfrac=mfrac,drift=drift,down=down,
                      group="sample"),
        random_sample(3,4;
                      nblk=nblk,nsig=nsig,sname="NIST612_1",
                      dtime=DateTime("2025-01-01T08:02:00"),
                      Time=Time,t0=t0,bwin=bwin,swin=swin,blank=blank,
                      mu_D=1000,sigma_D=50,y0=y0_glass,
                      p=0.0,mfrac=mfrac,group="NIST612"),
        random_sample(4,4;
                      nblk=nblk,nsig=nsig,sname="BP_2",
                      dtime=DateTime("2025-01-01T08:03:00"),
                      Time=Time,t0=t0,bwin=bwin,swin=swin,blank=blank,
                      mu_D=1000,sigma_D=50,x0=x0_std,y0=y0_std,
                      p=0.5,mfrac=mfrac,drift=drift,down=down,
                      group="BP_gt")
    ]
    return run
end
