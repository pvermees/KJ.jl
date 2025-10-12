# isochron
function getP(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
              vP::Real,vD::Real,vd::Real,
              sPD::Real,sPd::Real,sDd::Real,
              x0::Real,y0::Real,y1::Real,
              ft::AbstractVector,FT::AbstractVector,mf::Real,
              bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    return @. ((((bDt-Dm)*mf*vP+(Pm-bPt)*mf*sPD)*x0*y0+((dm-bdt)*mf^2*vP+(bPt-Pm)*mf^2*sPd)*x0)*y1+(((FT*Pm-FT*bPt)*ft*vD+(FT*bDt-Dm*FT)*ft*sPD)*x0^2+((Dm-bDt)*mf*vP+(bPt-Pm)*mf*sPD)*x0)*y0^2+(((Dm*FT-FT*bDt)*ft*mf*sPd+(FT*dm-FT*bdt)*ft*mf*sPD+(2*FT*bPt-2*FT*Pm)*ft*mf*sDd)*x0^2+((bdt-dm)*mf^2*vP+(Pm-bPt)*mf^2*sPd)*x0)*y0+((FT*Pm-FT*bPt)*ft*mf^2*vd+(FT*bdt-FT*dm)*ft*mf^2*sPd)*x0^2)/(mf^2*vP*y1^2+((2*FT*ft*mf*sPD*x0-2*mf^2*vP)*y0-2*FT*ft*mf^2*sPd*x0)*y1+(FT^2*ft^2*vD*x0^2-2*FT*ft*mf*sPD*x0+mf^2*vP)*y0^2+(2*FT*ft*mf^2*sPd*x0-2*FT^2*ft^2*mf*sDd*x0^2)*y0+FT^2*ft^2*mf^2*vd*x0^2)
end
export getP
# isochron
function getD(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
              vP::Real,vD::Real,vd::Real,
              sPD::Real,sPd::Real,sDd::Real,
              x0::Real,y0::Real,y1::Real,
              ft::AbstractVector,FT::AbstractVector,mf::Real,
              bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    return @. -((((bDt-Dm)*mf*vP+(Pm-bPt)*mf*sPD)*y1^2+((((FT*Pm-FT*bPt)*ft*vD+(FT*bDt-Dm*FT)*ft*sPD)*x0+(2*Dm-2*bDt)*mf*vP+(2*bPt-2*Pm)*mf*sPD)*y0+((2*Dm*FT-2*FT*bDt)*ft*mf*sPd+(FT*bdt-FT*dm)*ft*mf*sPD+(FT*bPt-FT*Pm)*ft*mf*sDd)*x0)*y1+(((FT*bPt-FT*Pm)*ft*vD+(Dm*FT-FT*bDt)*ft*sPD)*x0+(bDt-Dm)*mf*vP+(Pm-bPt)*mf*sPD)*y0^2+(((FT^2*bdt-FT^2*dm)*ft^2*vD+(Dm*FT^2-FT^2*bDt)*ft^2*sDd)*x0^2+((2*FT*bDt-2*Dm*FT)*ft*mf*sPd+(FT*dm-FT*bdt)*ft*mf*sPD+(FT*Pm-FT*bPt)*ft*mf*sDd)*x0)*y0+((FT^2*bDt-Dm*FT^2)*ft^2*mf*vd+(FT^2*dm-FT^2*bdt)*ft^2*mf*sDd)*x0^2)/(mf^2*vP*y1^2+((2*FT*ft*mf*sPD*x0-2*mf^2*vP)*y0-2*FT*ft*mf^2*sPd*x0)*y1+(FT^2*ft^2*vD*x0^2-2*FT*ft*mf*sPD*x0+mf^2*vP)*y0^2+(2*FT*ft*mf^2*sPd*x0-2*FT^2*ft^2*mf*sDd*x0^2)*y0+FT^2*ft^2*mf^2*vd*x0^2))
end
# point
function getD(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
              vP::Real,vD::Real,vd::Real,
              sPD::Real,sPd::Real,sDd::Real,
              x::Real,y::Real,
              ft::AbstractVector,FT::AbstractVector,mf::Real,
              bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    return @. ((((dm-bdt)*vD+(bDt-Dm)*sDd)*vP+(bPt-Pm)*sPd*vD+(Dm-bDt)*sPD*sPd+(bdt-dm)*sPD^2+(Pm-bPt)*sDd*sPD)*y+(((FT*Pm-FT*bPt)*ft*vD+(FT*bDt-Dm*FT)*ft*sPD)*vd+(FT*bdt-FT*dm)*ft*sPd*vD+(Dm*FT-FT*bDt)*ft*sDd*sPd+(FT*dm-FT*bdt)*ft*sDd*sPD+(FT*bPt-FT*Pm)*ft*sDd^2)*x+((Dm-bDt)*mf*vP+(bPt-Pm)*mf*sPD)*vd+(bdt-dm)*mf*sDd*vP+(bDt-Dm)*mf*sPd^2+((dm-bdt)*mf*sPD+(Pm-bPt)*mf*sDd)*sPd)/((vD*vP-sPD^2)*y^2+((2*FT*ft*sDd*sPD-2*FT*ft*sPd*vD)*x-2*mf*sDd*vP+2*mf*sPD*sPd)*y+(FT^2*ft^2*vD*vd-FT^2*ft^2*sDd^2)*x^2+(2*FT*ft*mf*sDd*sPd-2*FT*ft*mf*sPD*vd)*x+mf^2*vP*vd-mf^2*sPd^2)
end
# glass
function getD(Dm::AbstractVector,dm::AbstractVector,
              vD::Real,vd::Real,sDd::Real,
              y::Real,
              mf::Real,
              bDt::AbstractVector,bdt::AbstractVector)
    return (((dm-bdt)*vD+(bDt-Dm)*sDd)*y+(Dm-bDt)*mf*vd+(bdt-dm)*mf*sDd)/(vD*y^2-2*mf*sDd*y+mf^2*vd)
end
export getD

# mass fractionation + elemental fractionation
function SS(par::AbstractVector,
            bP::AbstractVector,bD::AbstractVector,bd::AbstractVector,
            dats::AbstractDict,covs::AbstractDict,
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
        covmat = covs[refmat]
        a = anchors[refmat]
        for spot in eachindex(dat)
            Pm,Dm,dm,vP,vD,vd,sPD,sPd,sDd,ft,FT,mf,bPt,bDt,bdt =
                SSprep(bP,bD,bd,dat[spot],covmat[spot],
                       channels,mfrac,drift,down;
                       PAcutoff=PAcutoff,adrift=adrift)
            if is_isochron_anchor(a)
                out += SS(Pm,Dm,dm,vP,vD,vd,sPD,sPd,sDd,
                          a.x0,a.y0,a.y1,ft,FT,mf,bPt,bDt,bdt)
            elseif is_point_anchor(a)
                out += SS(Pm,Dm,dm,vP,vD,vd,sPD,sPD,sDd,
                          a.x,a.y,ft,FT,mf,bPt,bDt,bdt)
            else
                error("Invalid anchor.")
            end
        end
    end
    return out
end
# isochron
function SS(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
            vP::Real,vD::Real,vd::Real,
            sPD::Real,sPd::Real,sDd::Real,
            x0::Real,y0::Real,y1::Real,
            ft::AbstractVector,FT::AbstractVector,mf::Real,
            bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    P = getP(Pm,Dm,dm,vP,vD,vd,sPD,sPd,sDd,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    D = getD(Pm,Dm,dm,vP,vD,vd,sPD,sPd,sDd,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    maha = @. (-((P*(y1-y0))/x0)-D*y0+dm-bdt)*(((vD*vP-sPD^2)*(-((P*(y1-y0))/x0)-D*y0+dm-bdt))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((-(D*mf)-bDt+Dm)*(sPD*sPd-sDd*vP))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((-(FT*P*ft)-bPt+Pm)*(sDd*sPD-sPd*vD))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD)))+(-(D*mf)-bDt+Dm)*(((sPD*sPd-sDd*vP)*(-((P*(y1-y0))/x0)-D*y0+dm-bdt))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((-(D*mf)-bDt+Dm)*(vP*vd-sPd^2))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((-(FT*P*ft)-bPt+Pm)*(sDd*sPd-sPD*vd))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD)))+(-(FT*P*ft)-bPt+Pm)*(((sDd*sPD-sPd*vD)*(-((P*(y1-y0))/x0)-D*y0+dm-bdt))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((-(FT*P*ft)-bPt+Pm)*(vD*vd-sDd^2))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((-(D*mf)-bDt+Dm)*(sDd*sPd-sPD*vd))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD)))
    return sum(@. maha )
end
# point
function SS(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
            vP::Real,vD::Real,vd::Real,
            sPD::Real,sPd::Real,sDd::Real,
            x::Real,y::Real,
            ft::AbstractVector,FT::AbstractVector,mf::Real,
            bPt::AbstractVector,bDt::AbstractVector,bdt::AbstractVector)
    D = getD(Pm,Dm,dm,vP,vD,vd,sPD,sPd,sDd,x,y,ft,FT,mf,bPt,bDt,bdt)
    maha = @. (-(D*y)+dm-bdt)*(((vD*vP-sPD^2)*(-(D*y)+dm-bdt))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((sDd*sPD-sPd*vD)*(-(D*FT*ft*x)-bPt+Pm))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((-(D*mf)-bDt+Dm)*(sPD*sPd-sDd*vP))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD)))+(-(D*mf)-bDt+Dm)*(((sPD*sPd-sDd*vP)*(-(D*y)+dm-bdt))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((sDd*sPd-sPD*vd)*(-(D*FT*ft*x)-bPt+Pm))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((-(D*mf)-bDt+Dm)*(vP*vd-sPd^2))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD)))+(-(D*FT*ft*x)-bPt+Pm)*(((sDd*sPD-sPd*vD)*(-(D*y)+dm-bdt))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((vD*vd-sDd^2)*(-(D*FT*ft*x)-bPt+Pm))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD))+((-(D*mf)-bDt+Dm)*(sDd*sPd-sPD*vd))/(vP*(vD*vd-sDd^2)+sPD*(sDd*sPd-sPD*vd)+sPd*(sDd*sPD-sPd*vD)))
    return sum(@. maha )
end
# glass
function SS(par::AbstractVector,
            bD::AbstractVector,bd::AbstractVector,
            dats::AbstractDict,covs::AbstractDict,
            channels::AbstractDict,
            anchors::AbstractDict)
    mf = exp(par[1])
    out = 0.0
    for (refmat,dat) in dats
        y = anchors[refmat]
        covmat = covs[refmat]
        for spot in eachindex(dat)
            Dm,dm,vD,vd,sDd,bDt,bdt =
                SSprep(bD,bd,dat[spot],covmat[spot],channels)
            out += SS(Dm,dm,vD,vd,sDd,y,mf,bDt,bdt)
        end
    end
    return out
end
# glass
function SS(Dm::AbstractVector,dm::AbstractVector,
            vD::Real,vd::Real,sDd::Real,
            y::Real,
            mf::Real,
            bDt::AbstractVector,bdt::AbstractVector)
    D = getD(Dm,dm,vD,vd,sDd,y,mf,bDt,bdt)
    maha = @. (-(D*y)+dm-bdt)*((vD*(-(D*y)+dm-bdt))/(vD*vd-sDd^2)-((-(D*mf)-bDt+Dm)*sDd)/(vD*vd-sDd^2))+(-(D*mf)-bDt+Dm)*(((-(D*mf)-bDt+Dm)*vd)/(vD*vd-sDd^2)-(sDd*(-(D*y)+dm-bdt))/(vD*vd-sDd^2))
    return sum(@. maha )
end
export SS

# isochron or point
function SSprep(bP::AbstractVector,bD::AbstractVector,bd::AbstractVector,
                dat::AbstractDataFrame,covmat::Matrix,
                channels::AbstractDict,mfrac::Real,
                drift::AbstractVector,down::AbstractVector;
                PAcutoff::Union{Real,Nothing}=nothing,
                adrift::AbstractVector=drift)
    t = dat.t
    T = dat.T
    sig = getSignals(dat)
    iP = columnindex(sig,channels["P"])
    iD = columnindex(sig,channels["D"])
    id = columnindex(sig,channels["d"])
    Pm = sig[:,iP]
    Dm = sig[:,iD]
    dm = sig[:,id]
    vP = covmat[iP,iP]
    vD = covmat[iD,iD]
    vd = covmat[id,id]
    sPD = covmat[iP,iD]
    sPd = covmat[iP,id]
    sDd = covmat[iD,id]
    ft = get_drift(Pm,t,drift;
                   PAcutoff=PAcutoff,adrift=adrift)
    FT = polyFac(down,T)
    mf = exp(mfrac)
    bPt = polyVal(bP,t)
    bDt = polyVal(bD,t)
    bdt = polyVal(bd,t)
    return Pm,Dm,dm,vP,vD,vd,sPD,sPd,sDd,ft,FT,mf,bPt,bDt,bdt
end
# glass
function SSprep(bD::AbstractVector,bd::AbstractVector,
                dat::AbstractDataFrame,covmat::Matrix,
                channels::AbstractDict)
    t = dat.t
    sig = getSignals(dat)
    iD = columnindex(sig,channels["D"])
    id = columnindex(sig,channels["d"])
    Dm = sig[:,iD]
    dm = sig[:,id]
    vD = covmat[iD,iD]
    vd = covmat[id,id]
    sDd = covmat[iD,id]
    bDt = polyVal(bD,t)
    bdt = polyVal(bd,t)
    return Dm,dm,vD,vd,sDd,bDt,bdt
end
export SSprep

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
        sig = getSignals(dat)
        covmat = df2cov(sig)
        anchor = anchors[samp.group]
        return predict(dat,covmat,pars,blank,channels,anchor;debug=debug)
    end
end
function predict(dat::AbstractDataFrame,
                 covmat::Matrix,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 anchor::NamedTuple;
                 debug::Bool=false)
    bP = blank[:,channels["P"]]
    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    Pm,Dm,dm,vP,vD,vd,sPD,sPd,sDd,ft,FT,mf,bPt,bDt,bdt =
        SSprep(bP,bD,bd,dat,covmat,channels,
               pars.mfrac,pars.drift,pars.down;
               PAcutoff=pars.PAcutoff,adrift=pars.adrift)
    if is_isochron_anchor(anchor)
        return predict(Pm,Dm,dm,vP,vD,vd,sPD,sPd,sDd,
                       anchor.x0,anchor.y0,anchor.y1,
                       ft,FT,mf,bPt,bDt,bdt)
    elseif is_point_anchor(anchor)
        return predict(Pm,Dm,dm,vP,vD,vd,sPD,sPd,sDd,
                       anchor.x,anchor.y,
                       ft,FT,mf,bPt,bDt,bdt)        
    else
        error("Invalid anchor")
    end
end
# isochron
function predict(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
                 vP::Real,vD::Real,vd::Real,
                 sPD::Real,sPd::Real,sDd::Real,
                 x0::Real,y0::Real,y1::Real,
                 ft::AbstractVector,FT::AbstractVector,mf::Real,
                 bPt::AbstractVector,bDt::AbstractVector,
                 bdt::AbstractVector)
    P = getP(Pm,Dm,dm,vP,vD,vd,sPD,sPd,sDd,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    D = getD(Pm,Dm,dm,vP,vD,vd,sPD,sPd,sDd,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    Pf = @. P*ft*FT + bPt
    Df = @. D*mf + bDt
    df = @. D*y0 + P*(y1-y0)/x0 + bdt
    return DataFrame(P=Pf,D=Df,d=df)
end
# point
function predict(Pm::AbstractVector,Dm::AbstractVector,dm::AbstractVector,
                 vP::Real,vD::Real,vd::Real,
                 sPD::Real,sPd::Real,sDd::Real,
                 x::Real,y::Real,
                 ft::AbstractVector,FT::AbstractVector,mf::Real,
                 bPt::AbstractVector,bDt::AbstractVector,
                 bdt::AbstractVector)
    D = getD(Pm,Dm,dm,vP,vD,vd,sPD,sPd,sDd,x,y,ft,FT,mf,bPt,bDt,bdt)
    Pf = @. D*x*ft*FT + bPt
    Df = @. D*mf + bDt
    df = @. D*y + bdt
    return DataFrame(P=Pf,D=Df,d=df)
end
# glass
function predict(dat::AbstractDataFrame,
                 covmat::Matrix,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 y::Real;
                 debug::Bool=false)
    mf = exp(pars.mfrac)
    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    Dm,dm,vD,vd,sDd,bDt,bdt = SSprep(bD,bd,dat,covmat,channels)
    return predict(Dm,dm,vD,vd,sDd,y,mf,bDt,bdt)
end
function predict(Dm::AbstractVector,dm::AbstractVector,
                 vD::Real,vd::Real,sDd::Real,
                 y::Real,
                 mf::Real,
                 bDt::AbstractVector,
                 bdt::AbstractVector)
    D = getD(Dm,dm,vD,vd,sDd,y,mf,bDt,bdt)
    Df = @. D*mf + bDt
    df = @. D*y + bdt
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
