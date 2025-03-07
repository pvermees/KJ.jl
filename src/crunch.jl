# minerals
function getP(Pm::AbstractVector,
              Dm::AbstractVector,
              dm::AbstractVector,
              vP::AbstractVector,
              vD::AbstractVector,
              vd::AbstractVector,
              x0::AbstractFloat,
              y0::AbstractFloat,
              y1::AbstractFloat,
              ft::AbstractVector,
              FT::AbstractVector,
              mf::AbstractFloat,
              bPt::AbstractVector,
              bDt::AbstractVector,
              bdt::AbstractVector)
    return @. -((Dm*FT*ft*mf*vP*x0*y0-FT*dm*ft*vP*x0)*y1+((-Pm*mf^2*vD*x0^2)-Dm*FT*ft*mf*vP*x0)*y0^2+FT*dm*ft*vP*x0*y0-Pm*vd*x0^2)/(FT^2*ft^2*vP*y1^2-2*FT^2*ft^2*vP*y0*y1+(mf^2*vD*x0^2+FT^2*ft^2*vP)*y0^2+vd*x0^2)
end
function getD(Pm::AbstractVector,
              Dm::AbstractVector,
              dm::AbstractVector,
              vP::AbstractVector,
              vD::AbstractVector,
              vd::AbstractVector,
              x0::AbstractFloat,
              y0::AbstractFloat,
              y1::AbstractFloat,
              ft::AbstractVector,
              FT::AbstractVector,
              mf::AbstractFloat,
              bPt::AbstractVector,
              bDt::AbstractVector,
              bdt::AbstractVector)
    return @. (Dm*FT^2*ft^2*vP*y1^2+((-FT*Pm*ft*mf*vD*x0)-2*Dm*FT^2*ft^2*vP)*y0*y1+(FT*Pm*ft*mf*vD*x0+Dm*FT^2*ft^2*vP)*y0^2+dm*mf*vD*x0^2*y0+Dm*vd*x0^2)/(FT^2*ft^2*vP*y1^2-2*FT^2*ft^2*vP*y0*y1+(mf^2*vD*x0^2+FT^2*ft^2*vP)*y0^2+vd*x0^2)
end
# glass
function getD(Dm::AbstractVector,
              dm::AbstractVector,
              vD::AbstractVector,
              vd::AbstractVector,
              y0::AbstractFloat,
              mf::AbstractFloat,
              bDt::AbstractVector,
              bdt::AbstractVector)
    return @. (dm*mf*vD*y0+Dm*vd)/(mf^2*vD*y0^2+vd)
end

# mineral
function SS(par::AbstractVector,
            bP::AbstractVector,
            bD::AbstractVector,
            bd::AbstractVector,
            dats::AbstractDict,
            channels::AbstractDict,
            anchors::AbstractDict,
            mf::Union{AbstractFloat,Nothing};
            dt::Union{AbstractDict,Nothing}=nothing,
            dead::AbstractFloat=0.0,
            ndrift::Integer=1,
            ndown::Integer=0,
            PAcutoff::Union{AbstractFloat,Nothing}=nothing)
    drift = par[1:ndrift]
    down = vcat(0.0,par[ndrift+1:ndrift+ndown])
    mfrac = isnothing(mf) ? par[ndrift+ndown+1] : log(mf)
    adrift = isnothing(PAcutoff) ? drift : par[end-ndrift+1:end]
    out = 0.0
    for (refmat,dat) in dats
        (x0,y0,y1) = anchors[refmat]
        Pm,Dm,dm,vP,vD,vd,ft,FT,mf,bPt,bDt,bdt =
            SSprep(bP,bD,bd,dat,channels,mfrac,drift,down;
                   dt=dt,dead=dead,PAcutoff=PAcutoff,adrift=adrift)
        out += SS(Pm,Dm,dm,vP,vD,vd,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    end
    return out
end
function SS(Pm::AbstractVector,
            Dm::AbstractVector,
            dm::AbstractVector,
            vP::AbstractVector,
            vD::AbstractVector,
            vd::AbstractVector,
            x0::AbstractFloat,
            y0::AbstractFloat,
            y1::AbstractFloat,
            ft::AbstractVector,
            FT::AbstractVector,
            mf::AbstractFloat,
            bPt::AbstractVector,
            bDt::AbstractVector,
            bdt::AbstractVector)
    pred = predict(Pm,Dm,dm,vP,vD,vd,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    dP = @. pred[:,"P"] - Pm
    dD = @. pred[:,"D"] - Dm
    dd = @. pred[:,"d"] - dm
    return sum(@. (dP^2)/vP + (dD^2)/vD + (dd^2)/vd )
end
# glass
function SS(par::AbstractVector,
            bD::AbstractVector,
            bd::AbstractVector,
            dats::AbstractDict,
            channels::AbstractDict,
            anchors::AbstractDict;
            dt::Union{AbstractDict,Nothing}=nothing,
            dead::AbstractFloat=0.0)
    mf = exp(par[1])
    out = 0.0
    for (refmat,dat) in dats
        y0 = anchors[refmat]
        Dm,dm,vD,vd,bDt,bdt = SSprep(bD,bd,dat,channels;
                                     dt=dt,dead=dead)
        out += SS(Dm,dm,vD,vd,y0,mf,bDt,bdt)
    end
    return out
end
function SS(Dm::AbstractVector,
            dm::AbstractVector,
            vD::AbstractVector,
            vd::AbstractVector,
            y0::AbstractFloat,
            mf::AbstractFloat,
            bDt::AbstractVector,
            bdt::AbstractVector)
    pred = predict(Dm,dm,vD,vd,y0,mf,bDt,bdt)
    dD = @. pred[:,"D"] - Dm
    dd = @. pred[:,"d"] - dm
    return sum(@. (dD^2)/vD + (dd^2)/vd )
end
export SS

# minerals
function SSprep(bP::AbstractVector,
                bD::AbstractVector,
                bd::AbstractVector,
                dat::AbstractDataFrame,
                channels::AbstractDict,
                mfrac::AbstractFloat,
                drift::AbstractVector,
                down::AbstractVector;
                dt::Union{AbstractDict,Nothing}=nothing,
                dead::AbstractFloat=0.0,
                PAcutoff::Union{AbstractFloat,Nothing}=nothing,
                adrift::AbstractVector=drift)
    t = dat.t
    T = dat.T
    Pm = dat[:,channels["P"]]
    Dm = dat[:,channels["D"]]
    dm = dat[:,channels["d"]]
    if isnothing(dt)
        vP = var_timeseries(Pm)
        vD = var_timeseries(Dm)
        vd = var_timeseries(dm)
    else
        vP = var_cps(Pm,dt[channels["P"]],dead)
        vD = var_cps(Dm,dt[channels["D"]],dead)
        vd = var_cps(dm,dt[channels["d"]],dead)
    end
    ft = get_drift(Pm,t,drift;
                   PAcutoff=PAcutoff,adrift=adrift)
    FT = polyFac(down,T)
    mf = exp(mfrac)
    bPt = polyVal(bP,t)
    bDt = polyVal(bD,t)
    bdt = polyVal(bd,t)
    return Pm,Dm,dm,vP,vD,vd,ft,FT,mf,bPt,bDt,bdt
end
# glass
function SSprep(bD::AbstractVector,
                bd::AbstractVector,
                dat::AbstractDataFrame,
                channels::AbstractDict;
                dt::Union{AbstractDict,Nothing}=nothing,
                dead::AbstractFloat=0.0)
    t = dat.t
    Dm = dat[:,channels["D"]]
    dm = dat[:,channels["d"]]
    if isnothing(dt)
        vD = var_timeseries(Dm)
        vd = var_timeseries(dm)
    else
        vD = var_cps(Dm,dt[channels["D"]],dead)
        vd = var_cps(dm,dt[channels["d"]],dead)
    end
    bDt = polyVal(bD,t)
    bdt = polyVal(bd,t)
    return Dm,dm,vD,vd,bDt,bdt
end

# minerals
function predict(samp::Sample,
                 method::AbstractString,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 standards::AbstractDict,
                 glass::AbstractDict;
                 dt::Union{AbstractDict,Nothing}=nothing,
                 dead::AbstractFloat=0.0,
                 debug::Bool=false)
    anchors = getAnchors(method,standards,glass)
    return predict(samp,pars,blank,channels,anchors;
                   dt=dt,dead=dead,debug=debug)
end
function predict(samp::Sample,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 anchors::AbstractDict;
                 dt::Union{AbstractDict,Nothing}=nothing,
                 dead::AbstractFloat=0.0,
                 debug::Bool=false)
    if samp.group == "sample"
        KJerror("notStandard")
    else
        dat = windowData(samp;signal=true)
        anchor = anchors[samp.group]
        return predict(dat,pars,blank,channels,anchor;
                       dt=dt,dead=dead,debug=debug)
    end
end
function predict(dat::AbstractDataFrame,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 anchor::NamedTuple;
                 dt::Union{AbstractDict,Nothing}=nothing,
                 dead::AbstractFloat=0.0,
                 debug::Bool=false)
    bP = blank[:,channels["P"]]
    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    Pm,Dm,dm,vP,vD,vd,ft,FT,mf,bPt,bDt,bdt =
        SSprep(bP,bD,bd,dat,channels,
               pars.mfrac,pars.drift,pars.down;
               dt=dt,dead=dead,
               PAcutoff=pars.PAcutoff,adrift=pars.adrift)
    return predict(Pm,Dm,dm,vP,vD,vd,
                   anchor.x0,anchor.y0,anchor.y1,
                   ft,FT,mf,bPt,bDt,bdt)
end
function predict(Pm::AbstractVector,
                 Dm::AbstractVector,
                 dm::AbstractVector,
                 vP::AbstractVector,
                 vD::AbstractVector,
                 vd::AbstractVector,
                 x0::AbstractFloat,
                 y0::AbstractFloat,
                 y1::AbstractFloat,
                 ft::AbstractVector,
                 FT::AbstractVector,
                 mf::AbstractFloat,
                 bPt::AbstractVector,
                 bDt::AbstractVector,
                 bdt::AbstractVector)
    P = getP(Pm,Dm,dm,vP,vD,vd,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    D = getD(Pm,Dm,dm,vP,vD,vd,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    return predict(P,D,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
end
function predict(P::AbstractVector,
                 D::AbstractVector,
                 x0::AbstractFloat,
                 y0::AbstractFloat,
                 y1::AbstractFloat,
                 ft::AbstractVector,
                 FT::AbstractVector,
                 mf::AbstractFloat,
                 bPt::AbstractVector,
                 bDt::AbstractVector,
                 bdt::AbstractVector)
    Pf = @. P + bPt
    Df = @. D + bDt
    df = @. D*y0*mf + P*ft*FT*(y1-y0)/x0 + bdt
    return DataFrame(P=Pf,D=Df,d=df)
end
# glass
function predict(dat::AbstractDataFrame,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 y0::AbstractFloat;
                 dt::Union{AbstractDict,Nothing}=nothing,
                 dead::AbstractFloat=0.0,
                 debug::Bool=false)
    mf = exp(pars.mfrac)
    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    Dm,dm,vD,vd,bDt,bdt = SSprep(bD,bd,dat,channels;
                                 dt=dt,dead=dead)
    return predict(Dm,dm,vD,vd,y0,mf,bDt,bdt)
end
function predict(Dm::AbstractVector,
                 dm::AbstractVector,
                 vD::AbstractVector,
                 vd::AbstractVector,
                 y0::AbstractFloat,
                 mf::AbstractFloat,
                 bDt::AbstractVector,
                 bdt::AbstractVector)
    D = getD(Dm,dm,vD,vd,y0,mf,bDt,bdt)
    Df = @. D + bDt
    df = @. D*mf*y0 + bdt
    return DataFrame(D=Df,d=df)
end
# concentrations
function predict(samp::Sample,
                 ef::AbstractVector,
                 blank::AbstractDataFrame,
                 elements::AbstractDataFrame,
                 internal::AbstractString;
                 debug::Bool=false)
    if samp.group in collect(keys(_KJ["glass"]))
        dat = windowData(samp;signal=true)
        sig = getSignals(dat)
        Xm = sig[:,Not(internal)]
        Sm = sig[:,internal]
        concs = elements2concs(elements,samp.group)
        R = collect((concs[:,Not(internal)]./concs[:,internal])[1,:])
        bt = polyVal(blank,dat.t)
        bXt = bt[:,Not(internal)]
        bSt = bt[:,internal]
        S = Sm.-bSt
        out = copy(sig)
        out[!,Not(internal)] = @. (R*ef)'*S + bXt
        return out
    else
        KJerror("notStandard")
    end
end
# blank
function predict(samp::Sample,
                 blank::AbstractDataFrame;
                 debug::Bool=false)
    dat = windowData(samp;blank=true)
    return polyVal(blank,dat.t)
end
export predict

function get_drift(Pm::AbstractVector,
                   t::AbstractVector,
                   pars::NamedTuple)
    return get_drift(Pm,t,pars.drift;
                     PAcutoff=pars.PAcutoff,
                     adrift=pars.adrift)
end
function get_drift(Pm::AbstractVector,
                   t::AbstractVector,
                   drift::AbstractVector;
                   PAcutoff=nothing,adrift=drift)
    if isnothing(PAcutoff)
        ft = polyFac(drift,t)
    else
        analog = Pm .> PAcutoff
        if all(analog)
            ft = polyFac(adrift,t)
        elseif all(.!analog)
            ft = polyFac(drift,t)
        else
            ft = polyFac(drift,t)
            ft[analog] = polyFac(adrift,t)[analog]
        end
    end
    return ft
end
