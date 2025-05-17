# isochron
function getP(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
              vP::AbstractVector,vD::AbstractVector,vd::AbstractVector,
              x0::AbstractFloat,y0::AbstractFloat,y1::AbstractFloat,
              ft::AbstractVector,FT::AbstractVector,mf::AbstractFloat,
              bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    return @. -((Dm*FT*ft*mf*vP*x0*y0-FT*dm*ft*vP*x0)*y1+((-Pm*mf^2*vD*x0^2)-Dm*FT*ft*mf*vP*x0)*y0^2+FT*dm*ft*vP*x0*y0-Pm*vd*x0^2)/(FT^2*ft^2*vP*y1^2-2*FT^2*ft^2*vP*y0*y1+(mf^2*vD*x0^2+FT^2*ft^2*vP)*y0^2+vd*x0^2)
end
# isochron
function getD(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
              vP::AbstractVector,vD::AbstractVector,vd::AbstractVector,
              x0::AbstractFloat,y0::AbstractFloat,y1::AbstractFloat,
              ft::AbstractVector,FT::AbstractVector,mf::AbstractFloat,
              bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    return @. (Dm*FT^2*ft^2*vP*y1^2+((-FT*Pm*ft*mf*vD*x0)-2*Dm*FT^2*ft^2*vP)*y0*y1+(FT*Pm*ft*mf*vD*x0+Dm*FT^2*ft^2*vP)*y0^2+dm*mf*vD*x0^2*y0+Dm*vd*x0^2)/(FT^2*ft^2*vP*y1^2-2*FT^2*ft^2*vP*y0*y1+(mf^2*vD*x0^2+FT^2*ft^2*vP)*y0^2+vd*x0^2)
end
# point
function getD(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
              vP::AbstractVector,vD::AbstractVector,vd::AbstractVector,
              x0::AbstractFloat,y0::AbstractFloat,
              ft::AbstractVector,FT::AbstractVector,mf::AbstractFloat,
              bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    return @. (dm*vD*vP*y0+(FT*Pm-FT*bPt)*ft*vD*vd*x0+Dm*mf*vP*vd)/(vD*vP*y0^2+FT^2*ft^2*vD*vd*x0^2+mf^2*vP*vd)
    return nothing
end
# glass
function getD(Dm::AbstractVector,dm::AbstractVector,
              vD::AbstractVector,vd::AbstractVector,
              y0::AbstractFloat,
              mf::AbstractFloat,
              bDt::AbstractVector,bdt::AbstractVector)
    return @. (dm*mf*vD*y0+Dm*vd)/(mf^2*vD*y0^2+vd)
end

# mass fractionation + elemental fractionation
function SS(par::AbstractVector,
            bP::AbstractVector,bD::AbstractVector,bd::AbstractVector,
            dats::AbstractDict,vars::AbstractDict,
            channels::AbstractDict,anchors::AbstractDict,
            mf::Union{AbstractFloat,Nothing};
            ndrift::Integer=1,
            ndown::Integer=0,
            PAcutoff::Union{AbstractFloat,Nothing}=nothing)
    drift = par[1:ndrift]
    down = vcat(0.0,par[ndrift+1:ndrift+ndown])
    mfrac = isnothing(mf) ? par[ndrift+ndown+1] : log(mf)
    adrift = isnothing(PAcutoff) ? drift : par[end-ndrift+1:end]
    out = 0.0
    for (refmat,dat) in dats
        Pm,Dm,dm,vP,vD,vd,ft,FT,mf,bPt,bDt,bdt =
            SSprep(bP,bD,bd,dat,vars[refmat],channels,mfrac,drift,down;
                   PAcutoff=PAcutoff,adrift=adrift)
        a = anchors[refmat]
        if a.type == "isochron"
            out += SS(Pm,Dm,dm,vP,vD,vd,a.x0,a.y0,a.y1,ft,FT,mf,bPt,bDt,bdt)
        else # "point"
            out += SS(Pm,Dm,dm,vP,vD,vd,a.x0,a.y0,ft,FT,mf,bPt,bDt,bdt)
        end
    end
    return out
end
# isochron
function SS(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
            vP::AbstractVector,vD::AbstractVector,vd::AbstractVector,
            x0::AbstractFloat,y0::AbstractFloat,y1::AbstractFloat,
            ft::AbstractVector,FT::AbstractVector,mf::AbstractFloat,
            bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    pred = predict(Pm,Dm,dm,vP,vD,vd,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    dP = @. pred[:,"P"] - Pm
    dD = @. pred[:,"D"] - Dm
    dd = @. pred[:,"d"] - dm
    return sum(@. (dP^2)/vP + (dD^2)/vD + (dd^2)/vd )
end
# point
function SS(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
            vP::AbstractVector,vD::AbstractVector,vd::AbstractVector,
            x0::AbstractFloat,y0::AbstractFloat,
            ft::AbstractVector,FT::AbstractVector,mf::AbstractFloat,
            bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    pred = predict(Pm,Dm,dm,vP,vD,vd,x0,y0,ft,FT,mf,bPt,bDt,bdt)
    dP = @. pred[:,"P"] - Pm
    dD = @. pred[:,"D"] - Dm
    dd = @. pred[:,"d"] - dm
    return sum(@. (dP^2)/vP + (dD^2)/vD + (dd^2)/vd )
end
# glass
function SS(par::AbstractVector,
            bD::AbstractVector,bd::AbstractVector,
            dats::AbstractDict,vars::AbstractDict,
            channels::AbstractDict,
            anchors::AbstractDict)
    mf = exp(par[1])
    out = 0.0
    for (refmat,dat) in dats
        y0 = anchors[refmat]
        Dm,dm,vD,vd,bDt,bdt = SSprep(bD,bd,dat,vars[refmat],channels)
        out += SS(Dm,dm,vD,vd,y0,mf,bDt,bdt)
    end
    return out
end
# glass
function SS(Dm::AbstractVector,dm::AbstractVector,
            vD::AbstractVector,vd::AbstractVector,
            y0::AbstractFloat,
            mf::AbstractFloat,
            bDt::AbstractVector,bdt::AbstractVector)
    pred = predict(Dm,dm,vD,vd,y0,mf,bDt,bdt)
    dD = @. pred[:,"D"] - Dm
    dd = @. pred[:,"d"] - dm
    return sum(@. (dD^2)/vD + (dd^2)/vd )
end
export SS

# isochron or point
function SSprep(bP::AbstractVector,bD::AbstractVector,bd::AbstractVector,
                dat::AbstractDataFrame,var::AbstractDataFrame,channels::AbstractDict,
                mfrac::AbstractFloat,drift::AbstractVector,down::AbstractVector;
                PAcutoff::Union{AbstractFloat,Nothing}=nothing,
                adrift::AbstractVector=drift)
    t = dat.t
    T = dat.T
    Pm = dat[:,channels["P"]]
    Dm = dat[:,channels["D"]]
    dm = dat[:,channels["d"]]
    vP = var[:,channels["P"]]
    vD = var[:,channels["D"]]
    vd = var[:,channels["d"]]
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
function SSprep(bD::AbstractVector,bd::AbstractVector,
                dat::AbstractDataFrame,var::AbstractDataFrame,
                channels::AbstractDict)
    t = dat.t
    Dm = dat[:,channels["D"]]
    dm = dat[:,channels["d"]]
    vD = var[:,channels["D"]]
    vd = var[:,channels["d"]]
    bDt = polyVal(bD,t)
    bdt = polyVal(bd,t)
    return Dm,dm,vD,vd,bDt,bdt
end

# isochron or point
function predict(samp::Sample,
                 method::AbstractString,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 standards::AbstractDict,
                 glass::AbstractDict;
                 debug::Bool=false)
    anchors = getAnchors(method,standards,glass)
    return predict(samp,pars,blank,channels,anchors;debug=debug)
end
function predict(samp::Sample,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 anchors::AbstractDict;
                 debug::Bool=false)
    if samp.group == "sample"
        KJerror("notStandard")
    else
        dat = windowData(samp;signal=true)
        var = dat2var(dat,collect(values(channels)))
        anchor = anchors[samp.group]
        return predict(dat,var,pars,blank,channels,anchor;debug=debug)
    end
end
function predict(dat::AbstractDataFrame,
                 var::AbstractDataFrame,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 anchor::NamedTuple;
                 debug::Bool=false)
    bP = blank[:,channels["P"]]
    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    Pm,Dm,dm,vP,vD,vd,ft,FT,mf,bPt,bDt,bdt =
        SSprep(bP,bD,bd,dat,var,channels,
               pars.mfrac,pars.drift,pars.down;
               PAcutoff=pars.PAcutoff,adrift=pars.adrift)
    if anchor.type == "isochron"
        return predict(Pm,Dm,dm,vP,vD,vd,
                       anchor.x0,anchor.y0,anchor.y1,
                       ft,FT,mf,bPt,bDt,bdt)
    else
        return predict(Pm,Dm,dm,vP,vD,vd,
                       anchor.x0,anchor.y0,
                       ft,FT,mf,bPt,bDt,bdt)        
    end
end
# isochron
function predict(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
                 vP::AbstractVector,vD::AbstractVector,vd::AbstractVector,
                 x0::AbstractFloat,y0::AbstractFloat,y1::AbstractFloat,
                 ft::AbstractVector,FT::AbstractVector,mf::AbstractFloat,
                 bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    P = getP(Pm,Dm,dm,vP,vD,vd,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    D = getD(Pm,Dm,dm,vP,vD,vd,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    Pf = @. P + bPt
    Df = @. D + bDt
    df = @. D*y0*mf + P*ft*FT*(y1-y0)/x0 + bdt
    return DataFrame(P=Pf,D=Df,d=df)
end
# point
function predict(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
                 vP::AbstractVector,vD::AbstractVector,vd::AbstractVector,
                 x0::AbstractFloat,y0::AbstractFloat,
                 ft::AbstractVector,FT::AbstractVector,mf::AbstractFloat,
                 bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    D = getD(Pm,Dm,dm,vP,vD,vd,x0,y0,ft,FT,mf,bPt,bDt,bdt)
    Pf = @. D*x0*ft*FT + bPt
    Df = @. D + bDt
    df = @. D*y0*mf + bdt
    return DataFrame(P=Pf,D=Df,d=df)
end
# glass
function predict(dat::AbstractDataFrame,
                 var::AbstractDataFrame,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 y0::AbstractFloat;
                 debug::Bool=false)
    mf = exp(pars.mfrac)
    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    Dm,dm,vD,vd,bDt,bdt = SSprep(bD,bd,dat,var,channels)
    return predict(Dm,dm,vD,vd,y0,mf,bDt,bdt)
end
function predict(Dm::AbstractVector,dm::AbstractVector,
                 vD::AbstractVector,vd::AbstractVector,
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
