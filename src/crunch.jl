function getP(a::IsochronAnchor,
              c::Cruncher,
              ft::AbstractVector,
              FT::AbstractVector)
    x0, y0, y1 = unpack(a)
    pmb, Dombi, bomb, bpt, bDot, bbot, vp, vD, vb, spD, spb, sDb, bd, t, T = unpack(c)
    return @. -((((Dombi*bd^2*vp-bd^2*pmb*spD)*x0*y0+(bd*pmb*spb-bd*bomb*vp)*x0)*y1+((Dombi*FT*bd^2*ft*spD-FT*bd^2*ft*pmb*vD)*x0^2+(bd^2*pmb*spD-Dombi*bd^2*vp)*x0)*y0^2+((-(Dombi*FT*bd*ft*spb)-FT*bd*bomb*ft*spD+2*FT*bd*ft*pmb*sDb)*x0^2+(bd*bomb*vp-bd*pmb*spb)*x0)*y0+(FT*bomb*ft*spb-FT*ft*pmb*vb)*x0^2)/(bd^2*vp*y1^2+((2*FT*bd^2*ft*spD*x0-2*bd^2*vp)*y0-2*FT*bd*ft*spb*x0)*y1+(FT^2*bd^2*ft^2*vD*x0^2-2*FT*bd^2*ft*spD*x0+bd^2*vp)*y0^2+(2*FT*bd*ft*spb*x0-2*FT^2*bd*ft^2*sDb*x0^2)*y0+FT^2*ft^2*vb*x0^2))
end
export getP

function getD(a::IsochronAnchor,
              c::Cruncher,
              ft::AbstractVector,
              FT::AbstractVector)
    x0, y0, y1 = unpack(a)
    pmb, Dombi, bomb, bpt, bDot, bbot, vp, vD, vb, spD, spb, sDb, bd, t, T = unpack(c)
    return @. ((Dombi*bd^2*vp-bd^2*pmb*spD)*y1^2+(((Dombi*FT*bd^2*ft*spD-FT*bd^2*ft*pmb*vD)*x0-2*Dombi*bd^2*vp+2*bd^2*pmb*spD)*y0+(-(2*Dombi*FT*bd*ft*spb)+FT*bd*bomb*ft*spD+FT*bd*ft*pmb*sDb)*x0)*y1+((FT*bd^2*ft*pmb*vD-Dombi*FT*bd^2*ft*spD)*x0+Dombi*bd^2*vp-bd^2*pmb*spD)*y0^2+((FT^2*bd*bomb*ft^2*vD-Dombi*FT^2*bd*ft^2*sDb)*x0^2+(2*Dombi*FT*bd*ft*spb-FT*bd*bomb*ft*spD-FT*bd*ft*pmb*sDb)*x0)*y0+(Dombi*FT^2*ft^2*vb-FT^2*bomb*ft^2*sDb)*x0^2)/(bd^2*vp*y1^2+((2*FT*bd^2*ft*spD*x0-2*bd^2*vp)*y0-2*FT*bd*ft*spb*x0)*y1+(FT^2*bd^2*ft^2*vD*x0^2-2*FT*bd^2*ft*spD*x0+bd^2*vp)*y0^2+(2*FT*bd*ft*spb*x0-2*FT^2*bd*ft^2*sDb*x0^2)*y0+FT^2*ft^2*vb*x0^2)
end

function getD(a::PointAnchor,
              c::Cruncher,
              ft::AbstractVector,
              FT::AbstractVector)
    x, y = unpack(a)
    pmb, Dombi, bomb, bpt, bDot, bbot, vp, vD, vb, spD, spb, sDb, bd, t, T = unpack(c)
    return @. (((bd*bomb*vD-Dombi*bd*sDb)*vp-bd*pmb*spb*vD+Dombi*bd*spD*spb-bd*bomb*spD^2+bd*pmb*sDb*spD)*y+((FT*ft*pmb*vD-Dombi*FT*ft*spD)*vb-FT*bomb*ft*spb*vD+Dombi*FT*ft*sDb*spb+FT*bomb*ft*sDb*spD-FT*ft*pmb*sDb^2)*x+(Dombi*vb-bomb*sDb)*vp-pmb*spD*vb-Dombi*spb^2+(bomb*spD+pmb*sDb)*spb)/((bd^2*vD*vp-bd^2*spD^2)*y^2+((2*FT*bd*ft*sDb*spD-2*FT*bd*ft*spb*vD)*x-2*bd*sDb*vp+2*bd*spD*spb)*y+(FT^2*ft^2*vD*vb-FT^2*ft^2*sDb^2)*x^2+(2*FT*ft*sDb*spb-2*FT*ft*spD*vb)*x+vb*vp-spb^2)
end

function SS(a::IsochronAnchor,
            c::Cruncher,
            ft::AbstractVector,
            FT::AbstractVector)
    Po = getP(a,c,ft,FT)
    Do = getD(a,c,ft,FT)
    x0, y0, y1 = unpack(a)
    pmb, Dombi, bomb, bpt, bDot, bbot, vp, vD, vb, spD, spb, sDb, bd, t, T = unpack(c)
    maha = @. (bomb-bd*((Po*(y1-y0))/x0+Do*y0))*(((vD*vp-spD^2)*(bomb-bd*((Po*(y1-y0))/x0+Do*y0)))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(spD*spb-sDb*vp))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((pmb-FT*Po*ft)*(sDb*spD-spb*vD))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(Dombi-Do)*(((spD*spb-sDb*vp)*(bomb-bd*((Po*(y1-y0))/x0+Do*y0)))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(vb*vp-spb^2))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((pmb-FT*Po*ft)*(sDb*spb-spD*vb))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(pmb-FT*Po*ft)*(((sDb*spD-spb*vD)*(bomb-bd*((Po*(y1-y0))/x0+Do*y0)))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((pmb-FT*Po*ft)*(vD*vb-sDb^2))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(sDb*spb-spD*vb))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))
    return sum(@. maha )
end

function SS(a::PointAnchor,
            c::Cruncher,
            ft::AbstractVector,
            FT::AbstractVector)
    Do = getD(a,c,ft,FT)
    x, y = unpack(a)
    pmb, Dombi, bomb, bpt, bDot, bbot, vp, vD, vb, spD, spb, sDb, bd, t, T = unpack(c)
    maha = @. (bomb-Do*bd*y)*(((vD*vp-spD^2)*(bomb-Do*bd*y))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((sDb*spD-spb*vD)*(pmb-Do*FT*ft*x))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(spD*spb-sDb*vp))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(Dombi-Do)*(((spD*spb-sDb*vp)*(bomb-Do*bd*y))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((sDb*spb-spD*vb)*(pmb-Do*FT*ft*x))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(vb*vp-spb^2))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(pmb-Do*FT*ft*x)*(((sDb*spD-spb*vD)*(bomb-Do*bd*y))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((vD*vb-sDb^2)*(pmb-Do*FT*ft*x))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(sDb*spb-spD*vb))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))
    
    return sum(@. maha )
    
end

function SS(par::AbstractVector,
            method::KJmethod,
            cruncher_groups::AbstractDict;
            verbose::Bool=false)

    fit = par2fit(par,method)
    out = 0.0
    for (refmat,crunchers) in cruncher_groups
        a = crunchers.anchor
        for c in crunchers.crunchers
            ft, FT = ft_FT(c,method,fit)
            out += SS(a,c,ft,FT)
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
        a = getAnchor(method.name,samp.group)
        c = Cruncher(samp,method,fit)
        ft, FT = ft_FT(c,method,fit)
        return predict(a,c,ft,FT)
    else
        KJerror("notStandard")
    end
end

function predict(a::IsochronAnchor,
                 c::Cruncher,
                 ft::AbstractVector,
                 FT::AbstractVector)
    Po = getP(a,c,ft,FT)
    Do = getD(a,c,ft,FT)
    x0, y0, y1 = unpack(a)
    pmb, Dombi, bomb, bpt, bDot, bbot, vp, vD, vb, spD, spb, sDb, bd, t, T = unpack(c)
    pf = @. Po*ft*FT + bpt
    Dof = @. Do + bDot
    bof = @. (Do*y0 + Po*(y1-y0)/x0)*bd + bbot
    return DataFrame(P=pf,D=Dof,d=bof)
end

function predict(a::PointAnchor,
                 c::Cruncher,
                 ft::AbstractVector,
                 FT::AbstractVector)
    Do = getD(a,c,ft,FT)
    x, y = unpack(a)
    pmb, Dombi, bomb, bpt, bDot, bbot, vp, vD, vb, spD, spb, sDb, bd, t, T = unpack(c)
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

function ft_FT(c::Cruncher,
               m::KJmethod,
               f::KJfit)
    ft = polyFac(f.drift,c.t)
    if !isnothing(m.PAcutoff)
        analog = c.pmb .> m.PAcutoff
        ft[analog] = polyFac(f.adrift,c.t)[analog]
    end
    FT = polyFac(f.down,c.T)
    return ft, FT
end

function Cruncher(samp::Sample,
                  method::KJmethod,
                  fit::KJfit)

    dat = swinData(samp)
    
    ch = getChannels(method)
    pm = dat[:,ch.P]
    Dom = dat[:,ch.D]
    bom = dat[:,ch.d]

    blank = fit.blank
    t = dat.t
    bpt = polyVal(blank[:,ch.P],t)
    bDot = polyVal(blank[:,ch.D],t)
    bbot = polyVal(blank[:,ch.d],t)

    pmb = pm - bpt
    Domb = Dom - bDot # TODO: add interference correction
    bomb = bom - bbot

    sig = hcat(pmb,Domb,bomb)
    covmat = df2cov(sig)
    vP = covmat[1,1]
    vD = covmat[2,2]
    vd = covmat[3,3]
    sPD = covmat[1,2]
    sPd = covmat[1,3]
    sDd = covmat[2,3]
    
    ions = getIons(method)
    proxies = getProxies(method)
    bd = iratio(proxies.d,ions.d)

    t = dat.t
    T = dat.T

    return Cruncher(pmb,Domb,bomb,bpt,bDot,bbot,vP,vD,vd,sPD,sPd,sDd,bd,t,T)
    
end
