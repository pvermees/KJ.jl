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
    P = getP(x0,y0,y1;
             pmb=pmb,Dombi=Dombi,bomb=bomb,
             vp=vp,vD=vD,vb=vb,spD=spD,spb=spb,sDb=sDb,
             ft=ft,FT=FT,bd=bd)
    D = getD(x0,y0,y1;
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
    D = getD(x,y;
             pmb=pmb,Dombi=Dombi,bomb=bomb,
             vp=vp,vD=vD,vb=vb,spD=spD,spb=spb,sDb=sDb,
             ft=ft,FT=FT,bd=bd)
    maha = @. (bomb-Do*bd*y)*(((vD*vp-spD^2)*(bomb-Do*bd*y))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((sDb*spD-spb*vD)*(pmb-Do*FT*ft*x))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(spD*spb-sDb*vp))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(Dombi-Do)*(((spD*spb-sDb*vp)*(bomb-Do*bd*y))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((sDb*spb-spD*vb)*(pmb-Do*FT*ft*x))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(vb*vp-spb^2))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(pmb-Do*FT*ft*x)*(((sDb*spD-spb*vD)*(bomb-Do*bd*y))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((vD*vb-sDb^2)*(pmb-Do*FT*ft*x))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(sDb*spb-spD*vb))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))
    return sum(@. maha )
end

"""
SS(par::NamedTuple,
   run::Vector{Sample},
   method::AbstractString,
   standards::AbstractDict,
   blank::AbstractDataFrame,
   channels::AbstractDict;
   verbose::Bool=false)

Sum of squares function for testing purposes
"""
function SS(par::NamedTuple,
            run::Vector{Sample},
            method::AbstractString,
            standards::AbstractDict,
            blank::AbstractDataFrame,
            channels::AbstractDict;
            verbose::Bool=false)
    anchors = getAnchors(method,standards)
    dats, covs, bP, bD, bd = SSfitprep(run,blank,anchors,channels)
    parvec = [par.drift; par.down]
    ndrift = length(par.drift)
    ndown = length(par.down)
    mf = exp(par.mfrac)
    return SS(parvec,bP,bD,bd,dats,covs,channels,anchors,mf;
              ndrift=ndrift,ndown=ndown,PAcutoff=par.PAcutoff,
              verbose=verbose)
end

"""
SS(par::AbstractVector,
   bP::AbstractVector,
   bD::AbstractVector,
   bd::AbstractVector,
   dats::AbstractDict,
   covs::AbstractDict,
   channels::AbstractDict,
   anchors::AbstractDict,
   mf::Union{Real,Nothing};
   ndrift::Integer=1,
   ndown::Integer=0,
   PAcutoff::Union{Real,Nothing}=nothing,
   verbose::Bool=false)

Sum of squares for mass fractionation + elemental fractionation
"""
function SS(par::AbstractVector,
            bP::AbstractVector,
            bD::AbstractVector,
            bd::AbstractVector,
            dats::AbstractDict,
            covs::AbstractDict,
            channels::AbstractDict,
            anchors::AbstractDict,
            mf::Union{Real,Nothing};
            ndrift::Integer=1,
            ndown::Integer=0,
            PAcutoff::Union{Real,Nothing}=nothing,
            verbose::Bool=false)
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
    if verbose
        println(par,": ",out)
    end
    return out
end
export SS

# isochron or point
function SSprep(dat::AbstractDataFrame,
                covmat::Matrix;
                bP::AbstractVector=[0.0],
                bD::AbstractVector=[0.0],
                bd::AbstractVector=[0.0],
                channels::AbstractDict,
                drift::AbstractVector,
                down::AbstractVector,
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

"""
predict(samp::Sample,
        method::KJmethod,
        fit::KJfit,
        anchors::AbstractDict)
"""
function predict(samp::Sample,
                 method::KJmethod,
                 fit::KJfit;
                 kw...)
    if samp.group in collect(keys(method.anchors))
        c = Cruncher(samp,method,fit)
        return predict(c)
    else
        KJerror("notStandard")
    end
end

function predict(c::Cruncher)
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
                       ft=c.ft,
                       FT=c.FT,
                       bd=c.bd)
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
                       ft=c.ft,
                       FT=c.FT,
                       bd=c.bd)
    else
        error("Invalid anchor")
    end
end

function predict(x0::Real,
                 y0::Real,
                 y1::Real=0.0;
                 pm::Vector{Real},
                 Dom::Vector{Real},
                 bom::Vector{Real},
                 bpt::Vector{Real},
                 bDot::Vector{Real},
                 bbot::Vector{Real},
                 vp::Real,
                 vD::Real,
                 vb::Real,
                 spD::Real,
                 spb::Real,
                 sDb::Real,
                 ft::Vector{Real},
                 FT::Vector{Real},
                 bd::Real)
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
              ft=ft,
              FT=FT,
              bd=bd)
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
              ft=ft,
              FT=FT,
              bd=bd)
    pf = @. Po*ft*FT + bpt
    Dof = @. Do + bDot
    bof = @. (Do*y0 + Po*(y1-y0)/x0)*bd + bbot
    return DataFrame(P=pf,D=Dof,d=bof)
end

function predict(x::Real,
                 y::Real;
                 pm::Vector{Real},
                 Dom::Vector{Real},
                 bom::Vector{Real},
                 bpt::Vector{Real},
                 bDot::Vector{Real},
                 bbot::Vector{Real},
                 vp::Real,
                 vD::Real,
                 vb::Real,
                 spD::Real,
                 spb::Real,
                 sDb::Real,
                 ft::Vector{Real},
                 FT::Vector{Real},
                 bd::Real)
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
              ft=ft,
              FT=FT,
              bd=bd)
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

    T = dat.T
    ft = get_drift(pm,t,fit)
    FT = polyFac(fit.down,T)
    
    ions = getIons(method)
    proxies = getProxies(method)
    bd = iratio(ions.d,proxies.d)

    return Cruncher(anchors,pm,Dom,bom,bpt,bDot,bbot,vP,vD,vd,sPD,sPd,sDd,ft,FT,bd)
end
