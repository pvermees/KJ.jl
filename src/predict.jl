function predict(samp::Sample,
                 method::Gmethod,
                 fit::Gfit)
    if samp.group in method.standards
        standard = method.groups[samp.group]
        a = getAnchor(method.name,standard)
        c = FCruncher(samp,method,fit)
        ft, hT = ft_hT(fit,method.PAcutoff;c...)
        return predict(a,ft,hT;c...)
    end
    for (element,calibration) in method.bias
        if samp.group in calibration.standards
            # TODO
        end
    end
    cg = bias_cruncher_groups_helper([samp],method.name,
                                     fit.blank,[samp.group],
                                     p.num,p.den)
    cruncher = cg[samp.group].crunchers[1]
    element = channel2element(p.num.channel)
    mf = polyFac(fit.bias[:,element],cruncher.t)
    y = cg[samp.group].anchor
    return predict(mf,y;cruncher...)
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
                 hT::AbstractVector,
                 P::AbstractVector,
                 D::AbstractVector;
                 bpt::AbstractVector,
                 bDt::AbstractVector,
                 bbt::AbstractVector,
                 mf::AbstractVector,
                 bd::AbstractFloat,
                 Ip::AbstractVector,
                 ID::AbstractVector,
                 Ib::AbstractVector,
                 other...)
    pf = @. P*ft*hT + bpt + Ip
    Df = @. D + bDt + ID
    bf = @. (D*a.y0 + P*(a.y1-a.y0)/a.x0)*bd*mf + bbt + Ib
    return DataFrame(P=pf,D=Df,d=bf)
end

function predict(a::PointAnchor,
                 ft::AbstractVector,
                 hT::AbstractVector,
                 D::AbstractVector;
                 bpt::AbstractVector,
                 bDt::AbstractVector,
                 bbt::AbstractVector,
                 bd::AbstractFloat,
                 mf::AbstractVector,
                 Ip::AbstractVector,
                 ID::AbstractVector,
                 Ib::AbstractVector,
                 other...)
    pf = @. D*a.x*ft*hT + bpt + Ip
    Df = @. D + bDt + ID
    bf = @. D*a.y*bd*mf + bbt + Ib
    return DataFrame(P=pf,D=Df,d=bf)
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
                 hT::AbstractVector;
                 cruncher...)
    P = getP(a,ft,hT;cruncher...)
    D = getD(a,ft,hT;cruncher...)
    return predict(a,ft,hT,P,D;cruncher...)
end

function predict(a::PointAnchor,
                 ft::AbstractVector,
                 hT::AbstractVector;
                 cruncher...)
    D = getD(a,ft,hT;cruncher...)
    return predict(a,ft,hT,D;cruncher...)
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

function ft_hT(f::Gfit,
               PAcutoff::AbstractFloat=Inf;
               pmb::AbstractVector,
               t::AbstractVector,
               T::AbstractVector,
               other...)
    ft = polyFac(f.drift,t)
    analog = pmb .> PAcutoff
    ft[analog] = polyFac(f.adrift,t)[analog]
    hT = polyFac(f.down,T)
    return ft, hT
end