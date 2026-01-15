function getP(a::IsochronAnchor,
              ft::AbstractVector,
              FT::AbstractVector;
              pmb::AbstractVector,
              Dombi::AbstractVector,
              bomb::AbstractVector,
              vp::AbstractFloat,
              vD::AbstractFloat,
              vb::AbstractFloat,
              spD::AbstractFloat,
              spb::AbstractFloat,
              sDb::AbstractFloat,
              bd::AbstractFloat,
              other...)
    x0,y0,y1 = unpack(a)
    return @. -((((Dombi*bd^2*vp-bd^2*pmb*spD)*x0*y0+(bd*pmb*spb-bd*bomb*vp)*x0)*y1+((Dombi*FT*bd^2*ft*spD-FT*bd^2*ft*pmb*vD)*x0^2+(bd^2*pmb*spD-Dombi*bd^2*vp)*x0)*y0^2+((-(Dombi*FT*bd*ft*spb)-FT*bd*bomb*ft*spD+2*FT*bd*ft*pmb*sDb)*x0^2+(bd*bomb*vp-bd*pmb*spb)*x0)*y0+(FT*bomb*ft*spb-FT*ft*pmb*vb)*x0^2)/(bd^2*vp*y1^2+((2*FT*bd^2*ft*spD*x0-2*bd^2*vp)*y0-2*FT*bd*ft*spb*x0)*y1+(FT^2*bd^2*ft^2*vD*x0^2-2*FT*bd^2*ft*spD*x0+bd^2*vp)*y0^2+(2*FT*bd*ft*spb*x0-2*FT^2*bd*ft^2*sDb*x0^2)*y0+FT^2*ft^2*vb*x0^2))
end
export getP

function getD(a::IsochronAnchor,
              ft::AbstractVector,
              FT::AbstractVector;
              pmb::AbstractVector,
              Dombi::AbstractVector,
              bomb::AbstractVector,
              vp::AbstractFloat,
              vD::AbstractFloat,
              vb::AbstractFloat,
              spD::AbstractFloat,
              spb::AbstractFloat,
              sDb::AbstractFloat,
              bd::AbstractFloat,
              other...)
    x0,y0,y1 = unpack(a)
    return @. ((Dombi*bd^2*vp-bd^2*pmb*spD)*y1^2+(((Dombi*FT*bd^2*ft*spD-FT*bd^2*ft*pmb*vD)*x0-2*Dombi*bd^2*vp+2*bd^2*pmb*spD)*y0+(-(2*Dombi*FT*bd*ft*spb)+FT*bd*bomb*ft*spD+FT*bd*ft*pmb*sDb)*x0)*y1+((FT*bd^2*ft*pmb*vD-Dombi*FT*bd^2*ft*spD)*x0+Dombi*bd^2*vp-bd^2*pmb*spD)*y0^2+((FT^2*bd*bomb*ft^2*vD-Dombi*FT^2*bd*ft^2*sDb)*x0^2+(2*Dombi*FT*bd*ft*spb-FT*bd*bomb*ft*spD-FT*bd*ft*pmb*sDb)*x0)*y0+(Dombi*FT^2*ft^2*vb-FT^2*bomb*ft^2*sDb)*x0^2)/(bd^2*vp*y1^2+((2*FT*bd^2*ft*spD*x0-2*bd^2*vp)*y0-2*FT*bd*ft*spb*x0)*y1+(FT^2*bd^2*ft^2*vD*x0^2-2*FT*bd^2*ft*spD*x0+bd^2*vp)*y0^2+(2*FT*bd*ft*spb*x0-2*FT^2*bd*ft^2*sDb*x0^2)*y0+FT^2*ft^2*vb*x0^2)
end

function getD(a::PointAnchor,
              ft::AbstractVector,
              FT::AbstractVector;
              pmb::AbstractVector,
              Dombi::AbstractVector,
              bomb::AbstractVector,
              vp::AbstractFloat,
              vD::AbstractFloat,
              vb::AbstractFloat,
              spD::AbstractFloat,
              spb::AbstractFloat,
              sDb::AbstractFloat,
              bd::AbstractFloat,
              other...)
    x,y = unpack(a)
    return @. (((bd*bomb*vD-Dombi*bd*sDb)*vp-bd*pmb*spb*vD+Dombi*bd*spD*spb-bd*bomb*spD^2+bd*pmb*sDb*spD)*y+((FT*ft*pmb*vD-Dombi*FT*ft*spD)*vb-FT*bomb*ft*spb*vD+Dombi*FT*ft*sDb*spb+FT*bomb*ft*sDb*spD-FT*ft*pmb*sDb^2)*x+(Dombi*vb-bomb*sDb)*vp-pmb*spD*vb-Dombi*spb^2+(bomb*spD+pmb*sDb)*spb)/((bd^2*vD*vp-bd^2*spD^2)*y^2+((2*FT*bd*ft*sDb*spD-2*FT*bd*ft*spb*vD)*x-2*bd*sDb*vp+2*bd*spD*spb)*y+(FT^2*ft^2*vD*vb-FT^2*ft^2*sDb^2)*x^2+(2*FT*ft*sDb*spb-2*FT*ft*spD*vb)*x+vb*vp-spb^2)
end

function getD(mf::AbstractVector,
              y::AbstractFloat;
              Dmb::AbstractVector,
              dmb::AbstractVector,
              vD::AbstractFloat,
              vd::AbstractFloat,
              sDd::AbstractFloat,
              other...)
    return @. ((dmb*mf*vD-Dmb*mf*sDd)*y+Dmb*vd-dmb*sDd)/(mf^2*vD*y^2-2*mf*sDd*y+vd)
end

function mahalanobis(a::IsochronAnchor,
                     ft::AbstractVector,
                     FT::AbstractVector,
                     Po::AbstractVector,
                     Do::AbstractVector;
                     pmb::AbstractVector,
                     Dombi::AbstractVector,
                     bomb::AbstractVector,
                     vp::AbstractFloat,
                     vD::AbstractFloat,
                     vb::AbstractFloat,
                     spD::AbstractFloat,
                     spb::AbstractFloat,
                     sDb::AbstractFloat,
                     bd::AbstractFloat,
                     other...)
    x0,y0,y1 = unpack(a)
    return @. (bomb-bd*((Po*(y1-y0))/x0+Do*y0))*(((vD*vp-spD^2)*(bomb-bd*((Po*(y1-y0))/x0+Do*y0)))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(spD*spb-sDb*vp))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((pmb-FT*Po*ft)*(sDb*spD-spb*vD))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(Dombi-Do)*(((spD*spb-sDb*vp)*(bomb-bd*((Po*(y1-y0))/x0+Do*y0)))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(vb*vp-spb^2))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((pmb-FT*Po*ft)*(sDb*spb-spD*vb))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(pmb-FT*Po*ft)*(((sDb*spD-spb*vD)*(bomb-bd*((Po*(y1-y0))/x0+Do*y0)))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((pmb-FT*Po*ft)*(vD*vb-sDb^2))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(sDb*spb-spD*vb))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))
end

function mahalanobis(a::PointAnchor,
                     ft::AbstractVector,
                     FT::AbstractVector,
                     Do::AbstractVector;
                     pmb::AbstractVector,
                     Dombi::AbstractVector,
                     bomb::AbstractVector,
                     vp::AbstractFloat,
                     vD::AbstractFloat,
                     vb::AbstractFloat,
                     spD::AbstractFloat,
                     spb::AbstractFloat,
                     sDb::AbstractFloat,
                     bd::AbstractFloat,
                     other...)
    x,y= unpack(a)
    return @. (bomb-Do*bd*y)*(((vD*vp-spD^2)*(bomb-Do*bd*y))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((sDb*spD-spb*vD)*(pmb-Do*FT*ft*x))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(spD*spb-sDb*vp))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(Dombi-Do)*(((spD*spb-sDb*vp)*(bomb-Do*bd*y))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((sDb*spb-spD*vb)*(pmb-Do*FT*ft*x))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(vb*vp-spb^2))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(pmb-Do*FT*ft*x)*(((sDb*spD-spb*vD)*(bomb-Do*bd*y))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((vD*vb-sDb^2)*(pmb-Do*FT*ft*x))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((Dombi-Do)*(sDb*spb-spD*vb))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))
end

function mahalanobis(mf::AbstractVector,
                     y::AbstractFloat,
                     D::AbstractVector;
                     Dmb::AbstractVector,
                     dmb::AbstractVector,
                     vD::AbstractFloat,
                     vd::AbstractFloat,
                     sDd::AbstractFloat,
                     other...)
    return @. (dmb-D*mf*y)*((vD*(dmb-D*mf*y))/(vD*vd-sDd^2)-((Dmb-D)*sDd)/(vD*vd-sDd^2))+(Dmb-D)*(((Dmb-D)*vd)/(vD*vd-sDd^2)-(sDd*(dmb-D*mf*y))/(vD*vd-sDd^2))
end

function SS(a::IsochronAnchor,
            ft::AbstractVector,
            FT::AbstractVector;
            cruncher...)
    Po = getP(a,ft,FT;cruncher...)
    Do = getD(a,ft,FT;cruncher...)
    maha = mahalanobis(a,ft,FT,Po,Do;cruncher...)
    return sum(@. maha )
end

function SS(a::PointAnchor,
            ft::AbstractVector,
            FT::AbstractVector;
            cruncher...)
    Do = getD(a,ft,FT;cruncher...)
    maha = mahalanobis(a,ft,FT,Do;cruncher...)
    return sum(@. maha )
end

function SS(par::AbstractVector,
            method::Gmethod,
            cruncher_groups::AbstractDict;
            verbose::Bool=false)
    fit = par2fit(par,method)
    out = 0.0
    for crunchers in values(cruncher_groups)
        a = crunchers.anchor
        for cruncher in crunchers.crunchers
            ft, FT = ft_FT(fit,method.PAcutoff;cruncher...)
            out += SS(a,ft,FT;cruncher...)
        end
    end
    if verbose
        println(par,": ",out)
    end
    return out
end

function SS(mf::AbstractVector,
            y::AbstractFloat;
            cruncher...)
    D = getD(mf,y;cruncher...)
    maha = mahalanobis(mf,y,D;cruncher...)
    return sum(@. maha )
end

function SS(par::AbstractVector,
            cruncher_groups::AbstractDict;
            verbose::Bool=false)
    out = 0.0
    # loop through standards:
    for cruncher_group in values(cruncher_groups)
        y = cruncher_group.anchor
        crunchers = cruncher_group.crunchers
        for cruncher in crunchers
            mf = polyFac(par,cruncher.t)
            out += SS(mf,y;cruncher...)
        end
    end
    if verbose
        println(par,": ",out)
    end
    return out
end

export SS

function predict(samp::Sample,
                 method::Gmethod,
                 fit::Gfit)
    if samp.group in getStandards(method.fractionation)
        a = getAnchor(method.name,samp.group)
        c = Cruncher(samp,method.fractionation,fit.blank)
        ft, FT = ft_FT(fit,method.PAcutoff;c...)
        return predict(a,ft,FT;c...)
    elseif samp.group in getStandards(method.fractionation.bias)
        p = bias_prep(method.fractionation)
        cg = bias_cruncher_groups_helper([samp],method.name,
                                         fit.blank,[samp.group],
                                         p.num,p.den)
        cruncher = cg[samp.group].crunchers[1]
        element = channel2element(p.num.channel)
        mf = polyFac(fit.bias[:,element],cruncher.t)
        y = cg[samp.group].anchor
        return predict(mf,y;cruncher...)
    elseif samp.group in getStandards(method.interference.bias)
        p = bias_prep(method.interference,target_channel,interfering_ion)
        c = Cruncher(samp,p.num,p.den,fit.blank)
        @infiltrate
    else
        KJerror("notStandard")
    end
end

function predict(samp::Sample,
                 method::Cmethod,
                 fit::Cfit)
    if samp.group in method.standards
        internal = method.internal[1]
        dat = swinData(samp)
        bt = polyVal(fit.blank,dat.t)
        S = getSignals(dat)[:,internal] .- bt[:,internal]
        C = getConcentrations(method,samp.group)
        Cs = C[1,internal]
        return @. S * fit.par * C / Cs + bt
    else
        KJerror("notStandard")
    end
end

function predict(a::IsochronAnchor,
                 ft::AbstractVector,
                 FT::AbstractVector,
                 Po::AbstractVector,
                 Do::AbstractVector;
                 bpt::AbstractVector,
                 bDot::AbstractVector,
                 bbot::AbstractVector,
                 bd::AbstractFloat,
                 other...)
    pf = @. Po*ft*FT + bpt
    Dof = @. Do + bDot
    bof = @. (Do*a.y0 + Po*(a.y1-a.y0)/a.x0)*bd + bbot
    return DataFrame(P=pf,D=Dof,d=bof)
end

function predict(a::PointAnchor,
                 ft::AbstractVector,
                 FT::AbstractVector,
                 Do::AbstractVector;
                 bpt::AbstractVector,
                 bDot::AbstractVector,
                 bbot::AbstractVector,
                 bd::AbstractFloat,
                 other...)
    pf = @. Do*a.x*ft*FT + bpt
    Dof = @. Do + bDot
    bof = @. Do*a.y*bd + bbot
    return DataFrame(P=pf,D=Dof,d=bof)
end

function predict(mf::AbstractVector,
                 y::AbstractFloat,
                 D::AbstractVector;
                 bDt::AbstractVector,
                 bdt::AbstractVector,
                 other...)
    Df = @. D + bDt
    df = @. D*y*mf + bdt
    return DataFrame(D=Df,d=df)
end

function predict(a::IsochronAnchor,
                 ft::AbstractVector,
                 FT::AbstractVector;
                 cruncher...)
    Po = getP(a,ft,FT;cruncher...)
    Do = getD(a,ft,FT;cruncher...)
    return predict(a,ft,FT,Po,Do;cruncher...)
end

function predict(a::PointAnchor,
                 ft::AbstractVector,
                 FT::AbstractVector;
                 cruncher...)
    Do = getD(a,ft,FT;cruncher...)
    return predict(a,ft,FT,Do;cruncher...)
end

function predict(mf::AbstractVector,
                 y::AbstractFloat;
                 cruncher...)
    D = getD(mf,y;cruncher...)
    return predict(mf,y,D;cruncher...)
end

function predict(samp::Sample,
                 ef::AbstractVector,
                 blank::AbstractDataFrame,
                 elements::AbstractDataFrame,
                 internal::AbstractString)
    if samp.group in _KJ["glass"].names
        dat = windowData(samp;signal=true)
        sig = getSignals(dat)
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

function predict(samp::Sample,
                 blank::AbstractDataFrame)
    dat = bwinData(samp)
    return polyVal(blank,dat.t)
end
export predict

function ft_FT(f::Gfit,
               PAcutoff::Union{Nothing,AbstractFloat};
               pmb::AbstractVector,
               t::AbstractVector,
               T::AbstractVector,
               other...)
    ft = polyFac(f.drift,t)
    if !isnothing(PAcutoff)
        analog = pmb .> PAcutoff
        ft[analog] = polyFac(f.adrift,t)[analog]
    end
    FT = polyFac(f.down,T)
    return ft, FT
end