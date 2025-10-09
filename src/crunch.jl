# isochron
function getP(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
              vP::AbstractVector,vD::AbstractVector,vd::AbstractVector,
              sPD::AbstractVector,sPd::AbstractVector,sDd::AbstractVector,
              x0::Real,y0::Real,y1::Real,
              ft::AbstractVector,FT::AbstractVector,mf::Real,
              bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    return @. -((((bDt-Dm)*m*vP+(Pm-bPt)*m*sPD)*x0*y0+((dm-bdt)*m^2*vP+(bPt-Pm)*m^2*sPd)*x0)*y1+(((FT*bPt-FT*Pm)*ft*vD+(Dm*FT-FT*bDt)*ft*sPD)*x0^2+((Dm-bDt)*m*vP+(bPt-Pm)*m*sPD)*x0)*y0^2+(((FT*bDt-Dm*FT)*ft*m*sPd+(FT*bdt-FT*dm)*ft*m*sPD+(2*FT*Pm-2*FT*bPt)*ft*m*sDd)*x0^2+((bdt-dm)*m^2*vP+(Pm-bPt)*m^2*sPd)*x0)*y0+((FT*bPt-FT*Pm)*ft*m^2*vd+(FT*dm-FT*bdt)*ft*m^2*sPd)*x0^2)/(m^2*vP*y1^2+(((-2*FT*ft*m*sPD*x0)-2*m^2*vP)*y0+2*FT*ft*m^2*sPd*x0)*y1+(FT^2*ft^2*vD*x0^2+2*FT*ft*m*sPD*x0+m^2*vP)*y0^2+((-2*FT^2*ft^2*m*sDd*x0^2)-2*FT*ft*m^2*sPd*x0)*y0+FT^2*ft^2*m^2*vd*x0^2)
end
# isochron
function getD(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
              vP::AbstractVector,vD::AbstractVector,vd::AbstractVector,
              sPD::AbstractVector,sPd::AbstractVector,sDd::AbstractVector,
              x0::Real,y0::Real,y1::Real,
              ft::AbstractVector,FT::AbstractVector,mf::Real,
              bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    return @. -(((bDt-Dm)*m*vP+(Pm-bPt)*m*sPD)*y1^2+((((FT*bPt-FT*Pm)*ft*vD+(Dm*FT-FT*bDt)*ft*sPD)*x0+(2*Dm-2*bDt)*m*vP+(2*bPt-2*Pm)*m*sPD)*y0+((2*FT*bDt-2*Dm*FT)*ft*m*sPd+(FT*dm-FT*bdt)*ft*m*sPD+(FT*Pm-FT*bPt)*ft*m*sDd)*x0)*y1+(((FT*Pm-FT*bPt)*ft*vD+(FT*bDt-Dm*FT)*ft*sPD)*x0+(bDt-Dm)*m*vP+(Pm-bPt)*m*sPD)*y0^2+(((FT^2*bdt-FT^2*dm)*ft^2*vD+(Dm*FT^2-FT^2*bDt)*ft^2*sDd)*x0^2+((2*Dm*FT-2*FT*bDt)*ft*m*sPd+(FT*bdt-FT*dm)*ft*m*sPD+(FT*bPt-FT*Pm)*ft*m*sDd)*x0)*y0+((FT^2*bDt-Dm*FT^2)*ft^2*m*vd+(FT^2*dm-FT^2*bdt)*ft^2*m*sDd)*x0^2)/(m^2*vP*y1^2+(((-2*FT*ft*m*sPD*x0)-2*m^2*vP)*y0+2*FT*ft*m^2*sPd*x0)*y1+(FT^2*ft^2*vD*x0^2+2*FT*ft*m*sPD*x0+m^2*vP)*y0^2+((-2*FT^2*ft^2*m*sDd*x0^2)-2*FT*ft*m^2*sPd*x0)*y0+FT^2*ft^2*m^2*vd*x0^2)
end
# point
function getD(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
              vP::AbstractVector,vD::AbstractVector,vd::AbstractVector,
              sPD::AbstractVector,sPd::AbstractVector,sDd::AbstractVector,
              x::Real,y::Real,
              ft::AbstractVector,FT::AbstractVector,mf::Real,
              bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    return @. ((((dm-bdt)*vD+(bDt-Dm)*sDd)*vP+(bPt-Pm)*sPd*vD+(Dm-bDt)*sPD*sPd+(bdt-dm)*sPD^2+(Pm-bPt)*sDd*sPD)*y+(((FT*Pm-FT*bPt)*ft*vD+(FT*bDt-Dm*FT)*ft*sPD)*vd+(FT*bdt-FT*dm)*ft*sPd*vD+(Dm*FT-FT*bDt)*ft*sDd*sPd+(FT*dm-FT*bdt)*ft*sDd*sPD+(FT*bPt-FT*Pm)*ft*sDd^2)*x+((Dm-bDt)*m*vP+(bPt-Pm)*m*sPD)*vd+(bdt-dm)*m*sDd*vP+(bDt-Dm)*m*sPd^2+((dm-bdt)*m*sPD+(Pm-bPt)*m*sDd)*sPd)/((vD*vP-sPD^2)*y^2+((2*FT*ft*sDd*sPD-2*FT*ft*sPd*vD)*x-2*m*sDd*vP+2*m*sPD*sPd)*y+(FT^2*ft^2*vD*vd-FT^2*ft^2*sDd^2)*x^2+(2*FT*ft*m*sDd*sPd-2*FT*ft*m*sPD*vd)*x+m^2*vP*vd-m^2*sPd^2)
end
# glass
function getD(Dm::AbstractVector,dm::AbstractVector,
              vD::AbstractVector,vd::AbstractVector,sDd::AbstractVector,
              y::Real,
              mf::Real,
              bDt::AbstractVector,bdt::AbstractVector)
    return @. (((dm-bdt)*vD+(bDt-Dm)*sDd)*y+(Dm-bDt)*m*vd+(bdt-dm)*m*sDd)/(vD*y^2-2*m*sDd*y+m^2*vd)
end

# mass fractionation + elemental fractionation
function SS(par::AbstractVector,
            bP::AbstractVector,bD::AbstractVector,bd::AbstractVector,
            dats::AbstractDict,vars::AbstractDict,
            channels::AbstractDict,anchors::AbstractDict,
            mf::Union{Real,Nothing};
            ndrift::Integer=1,
            ndown::Integer=0,
            PAcutoff::Union{Real,Nothing}=nothing)
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
        if isochronAnchor(a)
            out += SS(Pm,Dm,dm,vP,vD,vd,a.x0,a.y0,a.y1,ft,FT,mf,bPt,bDt,bdt)
        elseif pointAnchor(a)
            out += SS(Pm,Dm,dm,vP,vD,vd,a.x0,a.y0,ft,FT,mf,bPt,bDt,bdt)
        else
            error("Invalid anchor.")
        end
    end
    return out
end
# isochron
function SS(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
            vP::AbstractVector,vD::AbstractVector,vd::AbstractVector,
            x0::Real,y0::Real,y1::Real,
            ft::AbstractVector,FT::AbstractVector,mf::Real,
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
            x0::Real,y0::Real,
            ft::AbstractVector,FT::AbstractVector,mf::Real,
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
            y0::Real,
            mf::Real,
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
                mfrac::Real,drift::AbstractVector,down::AbstractVector;
                PAcutoff::Union{Real,Nothing}=nothing,
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
    if isochronAnchor(anchor)
        return predict(Pm,Dm,dm,vP,vD,vd,
                       anchor.x0,anchor.y0,anchor.y1,
                       ft,FT,mf,bPt,bDt,bdt)
    elseif pointAnchor(anchor)
        return predict(Pm,Dm,dm,vP,vD,vd,
                       anchor.x0,anchor.y0,
                       ft,FT,mf,bPt,bDt,bdt)        
    else
        error("Invalid anchor")
    end
end
# isochron
function predict(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
                 vP::AbstractVector,vD::AbstractVector,vd::AbstractVector,
                 x0::Real,y0::Real,y1::Real,
                 ft::AbstractVector,FT::AbstractVector,mf::Real,
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
                 x0::Real,y0::Real,
                 ft::AbstractVector,FT::AbstractVector,mf::Real,
                 bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    D = getD(Pm,Dm,dm,vP,vD,vd,x0,y0,ft,FT,mf,bPt,bDt,bdt)
    Pf = @. x0*D*mf/(ft*FT) + bPt
    Df = @. D + bDt
    df = @. y0*D*mf + bdt
    return DataFrame(P=Pf,D=Df,d=df)
end
# glass
function predict(dat::AbstractDataFrame,
                 var::AbstractDataFrame,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 y0::Real;
                 debug::Bool=false)
    mf = exp(pars.mfrac)
    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    Dm,dm,vD,vd,bDt,bdt = SSprep(bD,bd,dat,var,channels)
    return predict(Dm,dm,vD,vd,y0,mf,bDt,bdt)
end
function predict(Dm::AbstractVector,dm::AbstractVector,
                 vD::AbstractVector,vd::AbstractVector,
                 y0::Real,
                 mf::Real,
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
    if samp.group in _KJ["glass"].names
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
