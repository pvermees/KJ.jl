"""
getP(x0::Real,
     y0::Real,
     y1::Real;
     pmb::AbstractVector,
     Dombi::AbstractVector,
     bomb::AbstractVector,
     vp::Real,
     vD::Real,
     vb::Real,
     spD::Real,
     spb::Real,
     sDb::Real,
     ft::AbstractVector,
     FT::AbstractVector,
     bd::Real=1.0)

Estimate the true parent intensity from the measurements for isochron-based standards
"""
function getP(x0::Real,
              y0::Real,
              y1::Real;
              pmb::AbstractVector,
              Dombi::AbstractVector,
              bomb::AbstractVector,
              vp::Real,
              vD::Real,
              vb::Real,
              spD::Real,
              spb::Real,
              sDb::Real,
              ft::AbstractVector,
              FT::AbstractVector,
              bd::Real=1.0)
    return @. -((((Dombi*bd^2*vp-bd^2*pmb*spD)*x0*y0+(bd*pmb*spb-bd*bomb*vp)*x0)*y1+((Dombi*FT*bd^2*ft*spD-FT*bd^2*ft*pmb*vD)*x0^2+(bd^2*pmb*spD-Dombi*bd^2*vp)*x0)*y0^2+((-(Dombi*FT*bd*ft*spb)-FT*bd*bomb*ft*spD+2*FT*bd*ft*pmb*sDb)*x0^2+(bd*bomb*vp-bd*pmb*spb)*x0)*y0+(FT*bomb*ft*spb-FT*ft*pmb*vb)*x0^2)/(bd^2*vp*y1^2+((2*FT*bd^2*ft*spD*x0-2*bd^2*vp)*y0-2*FT*bd*ft*spb*x0)*y1+(FT^2*bd^2*ft^2*vD*x0^2-2*FT*bd^2*ft*spD*x0+bd^2*vp)*y0^2+(2*FT*bd*ft*spb*x0-2*FT^2*bd*ft^2*sDb*x0^2)*y0+FT^2*ft^2*vb*x0^2))
end
export getP

"""
getD(x0::Real,
     y0::Real,
     y1::Real;
     pmb::AbstractVector,
     Dombi::AbstractVector,
     bomb::AbstractVector,
     vp::Real,
     vD::Real,
     vb::Real,
     spD::Real,
     spb::Real,
     sDb::Real,
     ft::AbstractVector,
     FT::AbstractVector,
     bd::Number=1.0)

Estimate the true daughter intensity from the measurements for isochron-based standards
"""
function getD(x0::Real,
              y0::Real,
              y1::Real;
              pmb::AbstractVector,
              Dombi::AbstractVector,
              bomb::AbstractVector,
              vp::Real,
              vD::Real,
              vb::Real,
              spD::Real,
              spb::Real,
              sDb::Real,
              ft::AbstractVector,
              FT::AbstractVector,
              bd::Number=1.0)
    return @. ((Dombi*bd^2*vp-bd^2*pmb*spD)*y1^2+(((Dombi*FT*bd^2*ft*spD-FT*bd^2*ft*pmb*vD)*x0-2*Dombi*bd^2*vp+2*bd^2*pmb*spD)*y0+(-(2*Dombi*FT*bd*ft*spb)+FT*bd*bomb*ft*spD+FT*bd*ft*pmb*sDb)*x0)*y1+((FT*bd^2*ft*pmb*vD-Dombi*FT*bd^2*ft*spD)*x0+Dombi*bd^2*vp-bd^2*pmb*spD)*y0^2+((FT^2*bd*bomb*ft^2*vD-Dombi*FT^2*bd*ft^2*sDb)*x0^2+(2*Dombi*FT*bd*ft*spb-FT*bd*bomb*ft*spD-FT*bd*ft*pmb*sDb)*x0)*y0+(Dombi*FT^2*ft^2*vb-FT^2*bomb*ft^2*sDb)*x0^2)/(bd^2*vp*y1^2+((2*FT*bd^2*ft*spD*x0-2*bd^2*vp)*y0-2*FT*bd*ft*spb*x0)*y1+(FT^2*bd^2*ft^2*vD*x0^2-2*FT*bd^2*ft*spD*x0+bd^2*vp)*y0^2+(2*FT*bd*ft*spb*x0-2*FT^2*bd*ft^2*sDb*x0^2)*y0+FT^2*ft^2*vb*x0^2)
end

"""
getD(x::Real,
     y::Real;
     pmb::AbstractVector,
     Dombi::AbstractVector,
     bomb::AbstractVector,
     vp::Real,
     vD::Real,
     vb::Real,
     spD::Real,
     spb::Real,
     sDb::Real,
     ft::AbstractVector,
     FT::AbstractVector,
     bd::Number=1.0)

Estimate the true daughter intensity from the measurements for homogeneous standards
"""
function getD(x::Real,
              y::Real;
              pmb::AbstractVector,
              Dombi::AbstractVector,
              bomb::AbstractVector,
              vp::Real,
              vD::Real,
              vb::Real,
              spD::Real,
              spb::Real,
              sDb::Real,
              ft::AbstractVector,
              FT::AbstractVector,
              bd::Number=1.0)
    return @. (((bd*bomb*vD-Dombi*bd*sDb)*vp-bd*pmb*spb*vD+Dombi*bd*spD*spb-bd*bomb*spD^2+bd*pmb*sDb*spD)*y+((FT*ft*pmb*vD-Dombi*FT*ft*spD)*vb-FT*bomb*ft*spb*vD+Dombi*FT*ft*sDb*spb+FT*bomb*ft*sDb*spD-FT*ft*pmb*sDb^2)*x+(Dombi*vb-bomb*sDb)*vp-pmb*spD*vb-Dombi*spb^2+(bomb*spD+pmb*sDb)*spb)/((bd^2*vD*vp-bd^2*spD^2)*y^2+((2*FT*bd*ft*sDb*spD-2*FT*bd*ft*spb*vD)*x-2*bd*sDb*vp+2*bd*spD*spb)*y+(FT^2*ft^2*vD*vb-FT^2*ft^2*sDb^2)*x^2+(2*FT*ft*sDb*spb-2*FT*ft*spD*vb)*x+vb*vp-spb^2)
end

"""
SS(x0::Real,
   y0::Real,
   y1::Real;
   pmb::AbstractVector,
   Dombi::AbstractVector,
   bomb::AbstractVector,
   vp::Real,
   vD::Real,
   vb::Real,
   spD::Real,
   spb::Real,
   sDb::Real,
   ft::AbstractVector,
   FT::AbstractVector,
   bd::Number=1.0)

Sum of squares for isochron-based standards
"""
function SS(x0::Real,
            y0::Real,
            y1::Real;
            pmb::AbstractVector,
            Dombi::AbstractVector,
            bomb::AbstractVector,
            vp::Real,
            vD::Real,
            vb::Real,
            spD::Real,
            spb::Real,
            sDb::Real,
            ft::AbstractVector,
            FT::AbstractVector,
            bd::Number=1.0)
    Po = getP(x0,y0,y1;
              pmb=pmb,Dombi=Dombi,bomb=bomb,
              vp=vp,vD=vD,vb=vb,spD=spD,spb=spb,sDb=sDb,
              ft=ft,FT=FT,bd=bd)
    Do = getD(x0,y0,y1;
              pmb=pmb,Dombi=Dombi,bomb=bomb,
              vp=vp,vD=vD,vb=vb,spD=spD,spb=spb,sDb=sDb,
              ft=ft,FT=FT,bd=bd)
    maha = @. (bomb-bd*((Po*(y1-y0))/x0+Do*y0))*(((vD*vp-spD^2)*(bomb-bd*((Po*(y1-y0))/x0+Do*y0)))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(spD*spb-sDb*vp))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((pmb-FT*Po*ft)*(sDb*spD-spb*vD))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(Dombi-Do)*(((spD*spb-sDb*vp)*(bomb-bd*((Po*(y1-y0))/x0+Do*y0)))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(vb*vp-spb^2))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((pmb-FT*Po*ft)*(sDb*spb-spD*vb))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(pmb-FT*Po*ft)*(((sDb*spD-spb*vD)*(bomb-bd*((Po*(y1-y0))/x0+Do*y0)))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((pmb-FT*Po*ft)*(vD*vb-sDb^2))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(sDb*spb-spD*vb))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))
    return sum(@. maha )
end

"""
SS(x::Real,
   y::Real;
   pmb::AbstractVector,
   Dombi::AbstractVector,
   bomb::AbstractVector,
   vp::Real,
   vD::Real,
   vb::Real,
   spD::Real,
   spb::Real,
   sDb::Real,
   ft::AbstractVector,
   FT::AbstractVector,
   bd::Number=1.0)

Sum of squares for homogeneous standards
"""
function SS(x::Real,
            y::Real;
            pmb::AbstractVector,
            Dombi::AbstractVector,
            bomb::AbstractVector,
            vp::Real,
            vD::Real,
            vb::Real,
            spD::Real,
            spb::Real,
            sDb::Real,
            ft::AbstractVector,
            FT::AbstractVector,
            bd::Number=1.0)
    
    Do = getD(x,y;
              pmb=pmb,Dombi=Dombi,bomb=bomb,
              vp=vp,vD=vD,vb=vb,spD=spD,spb=spb,sDb=sDb,
              ft=ft,FT=FT,bd=bd)
    
    maha = @. (bomb-Do*bd*y)*(((vD*vp-spD^2)*(bomb-Do*bd*y))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((sDb*spD-spb*vD)*(pmb-Do*FT*ft*x))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(spD*spb-sDb*vp))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(Dombi-Do)*(((spD*spb-sDb*vp)*(bomb-Do*bd*y))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((sDb*spb-spD*vb)*(pmb-Do*FT*ft*x))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(vb*vp-spb^2))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(pmb-Do*FT*ft*x)*(((sDb*spD-spb*vD)*(bomb-Do*bd*y))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((vD*vb-sDb^2)*(pmb-Do*FT*ft*x))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(sDb*spb-spD*vb))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))
    
    return sum(@. maha )
    
end

function SS(par::AbstractVector,
            method::KJmethod,
            cruncher_groups::AbstractDict;
            verbose::Bool=false)

    drift, down, adrift = parse_par(par,method)
    
    out = 0.0
    for (refmat,crunchers) in cruncher_groups
        for c in crunchers
            
            ft = get_drift(c.pm,
                           c.t;
                           drift = drift,
                           adrift = adrift,
                           PAcutoff = method.PAcutoff)
            FT = polyFac(down,
                         c.T)
            a = c.anchor
            
            if is_isochron_anchor(a)
                out += SS(a.x0,
                          a.y0,
                          a.y1;
                          pmb = c.pm-c.bpt,
                          Dombi = c.Dom-c.bDot,
                          bomb = c.bom-c.bbot,
                          vp = c.vp,
                          vD = c.vD,
                          vb = c.vb,
                          spD = c.spD,
                          spb = c.spb,
                          sDb = c.sDb,
                          bd = c.bd,
                          ft = ft,
                          FT = FT)
            elseif is_point_anchor(a)
                out += SS(a.x,
                          a.y;
                          pmb = c.pm-c.bpt,
                          Dombi = c.Dom-c.bDot,
                          bomb = c.bom-c.bbot,
                          vp = c.vp,
                          vD = c.vD,
                          vb = c.vb,
                          spD = c.spD,
                          spb = c.spb,
                          sDb = c.sDb,
                          bd = c.bd,
                          ft = ft,
                          FT = FT)
            else
                error("Invalid anchor.")
            end
        end
    end
    if verbose
        println(par,": ",out)
    end
    return out
end
export SS

"""
predict(samp::Sample,
        method::KJmethod,
        fit::KJfit)
"""
function predict(samp::Sample,
                 method::KJmethod,
                 fit::KJfit;
                 kw...)
    if samp.group in collect(keys(method.anchors))
        return predict(fit.drift,
                       fit.down,
                       Cruncher(samp,method,fit);
                       PAcutoff=method.PAcutoff,
                       adrift=fit.adrift)
    else
        KJerror("notStandard")
    end
end

function predict(drift::AbstractVector,
                 down::AbstractVector,
                 c::Cruncher;
                 PAcutoff::Union{Nothing,Real}=nothing,
                 adrift::AbstractVector=drift)
    
    ft = get_drift(c.pm,
                   c.t;
                   drift=drift,
                   adrift=adrift,
                   PAcutoff=PAcutoff)
    FT = polyFac(down,
                 c.T)

    if is_isochron_anchor(c.anchor)
        return predict(c.anchor.x0,
                       c.anchor.y0,
                       c.anchor.y1;
                       pm=c.pm,
                       Dom=c.Dom,
                       bom=c.bom,
                       bpt=c.bpt,
                       bDot=c.bDot,
                       bbot=c.bbot,
                       vp=c.vp,
                       vD=c.vD,
                       vb=c.vb,
                       spD=c.spD,
                       spb=c.spb,
                       sDb=c.sDb,
                       bd=c.bd,
                       ft=ft,
                       FT=FT)
    elseif is_point_anchor(anchor)
        return predict(c.anchor.x,
                       c.anchor.y;
                       pm=c.pm,
                       Dom=c.Dom,
                       bom=c.bom,
                       bpt=c.bpt,
                       bDot=c.bDot,
                       bbot=c.bbot,
                       vp=c.vp,
                       vD=c.vD,
                       vb=c.vb,
                       spD=c.spD,
                       spb=c.spb,
                       sDb=c.sDb,
                       bd=c.bd,
                       ft=ft,
                       FT=FT)
    else
        error("Invalid anchor")
    end
end

function predict(x0::Real,
                 y0::Real,
                 y1::Real;
                 pm::AbstractVector,
                 Dom::AbstractVector,
                 bom::AbstractVector,
                 bpt::AbstractVector,
                 bDot::AbstractVector,
                 bbot::AbstractVector,
                 vp::Real,
                 vD::Real,
                 vb::Real,
                 spD::Real,
                 spb::Real,
                 sDb::Real,
                 bd::Real,
                 ft::AbstractVector,
                 FT::AbstractVector)
    Po = getP(x0,
              y0,
              y1;
              pmb=pm-bpt,
              Dombi=Dom-bDot,
              bomb=bom-bbot,
              vp=vp,
              vD=vD,
              vb=vb,
              spD=spD,
              spb=spb,
              sDb=sDb,
              bd=bd,
              ft=ft,
              FT=FT)
    Do = getD(x0,
              y0,
              y1;
              pmb=pm-bpt,
              Dombi=Dom-bDot,
              bomb=bom-bbot,
              vp=vp,
              vD=vD,
              vb=vb,
              spD=spD,
              spb=spb,
              sDb=sDb,
              bd=bd,
              ft=ft,
              FT=FT)
    pf = @. Po*ft*FT + bpt
    Dof = @. Do + bDot
    bof = @. (Do*y0 + Po*(y1-y0)/x0)*bd + bbot
    return DataFrame(P=pf,D=Dof,d=bof)
end

function predict(x::Real,
                 y::Real;
                 pm::AbstractVector,
                 Dom::AbstractVector,
                 bom::AbstractVector,
                 bpt::AbstractVector,
                 bDot::AbstractVector,
                 bbot::AbstractVector,
                 vp::Real,
                 vD::Real,
                 vb::Real,
                 spD::Real,
                 spb::Real,
                 sDb::Real,
                 bd::Real,
                 ft::AbstractVector,
                 FT::AbstractVector)
    Do = getD(x,
              y;
              pmb=pm-bpt,
              Dombi=Dom-bDot,
              bomb=bom-bbot,
              vp=vp,
              vD=vD,
              vb=vb,
              spD=spD,
              spb=spb,
              sDb=sDb,
              bd=bd,
              ft=ft,
              FT=FT)
    pf = @. Do*x*ft*FT + bpt
    Dof = @. Do + bDot
    bof = @. Do*y*bd + bbot
    return DataFrame(P=pf,D=Dof,d=bof)
end

"""
predict(samp::Sample,
        ef::AbstractVector,
        blank::AbstractDataFrame,
        elements::AbstractDataFrame,
        internal::AbstractString;
        debug::Bool=false)

For concentrations
"""
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

"""
predict(samp::Sample,
        blank::AbstractDataFrame;
        debug::Bool=false)

For blanks
"""
function predict(samp::Sample,
                 blank::AbstractDataFrame;
                 debug::Bool=false)
    dat = bwinData(samp)
    return polyVal(blank,dat.t)
end
export predict

function Cruncher(samp::Sample,
                  method::KJmethod,
                  fit::KJfit)

    dat = swinData(samp)
    sig = getSignals(dat)
    covmat = df2cov(sig)
    anchors = method.anchors[samp.group]
    
    ch = getChannels(method)
    pm = dat[:,ch.P]
    Dom = dat[:,ch.D]
    bom = dat[:,ch.d]

    blank = fit.blank
    t = dat.t
    bpt = polyVal(blank[:,ch.P],t)
    bDot = polyVal(blank[:,ch.D],t)
    bbot = polyVal(blank[:,ch.d],t)

    sig = getSignals(dat)
    iP = columnindex(sig,ch.P)
    iD = columnindex(sig,ch.D)
    id = columnindex(sig,ch.d)
    vP = covmat[iP,iP]
    vD = covmat[iD,iD]
    vd = covmat[id,id]
    sPD = covmat[iP,iD]
    sPd = covmat[iP,id]
    sDd = covmat[iD,id]
    
    ions = getIons(method)
    proxies = getProxies(method)
    bd = iratio(ions.d,proxies.d)

    t = dat.t
    T = dat.T

    return Cruncher(anchors,pm,Dom,bom,bpt,bDot,bbot,vP,vD,vd,sPD,sPd,sDd,bd,t,T)
    
end
