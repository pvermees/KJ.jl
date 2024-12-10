# for age standards
function w2S(wP::AbstractFloat,
             wd::AbstractFloat,
             Pm::AbstractVector,
             Dm::AbstractVector,
             dm::AbstractVector,
             x0::AbstractFloat,
             y0::AbstractFloat,
             y1::AbstractFloat,
             ft::AbstractVector,
             FT::AbstractVector,
             mf::AbstractFloat,
             bPt::AbstractVector,
             bDt::AbstractVector,
             bdt::AbstractVector)
    return @. (((dm-bdt)*mf^2+(Dm-bDt)*mf)*wd*y1^2+((((FT*bPt*ft-FT*ft*Pm)*wP+(bDt-Dm)*mf)*wd*x0+((2*bdt-2*dm)*mf^2+(2*bDt-2*Dm)*mf)*wd)*y0+((dm-bdt)*mf^2*wd+(FT*ft*Pm-FT*bPt*ft)*mf^2*wP)*x0)*y1+((FT^2*ft^2*dm+FT*ft*Pm-FT^2*bdt*ft^2-FT*bPt*ft)*wP*wd*x0^2+((FT*ft*Pm-FT*bPt*ft)*wP+(Dm-bDt)*mf)*wd*x0+((dm-bdt)*mf^2+(Dm-bDt)*mf)*wd)*y0^2+(((FT^2*ft^2*dm-FT^2*bdt*ft^2)*wP*wd+(FT^2*ft^2*Dm-FT^2*bDt*ft^2)*mf*wP)*x0^2+((bdt-dm)*mf^2*wd+(FT*bPt*ft-FT*ft*Pm)*mf^2*wP)*x0)*y0+((FT*ft*Pm-FT*bPt*ft)*mf^2+(FT^2*ft^2*Dm-FT^2*bDt*ft^2)*mf)*wP*x0^2)/(mf^2*wd*y1^2-2*mf^2*wd*y0*y1+(FT^2*ft^2*wP*wd*x0^2+mf^2*wd)*y0^2+FT^2*ft^2*mf^2*wP*x0^2)
end
function w2p(wP::AbstractFloat,
             wd::AbstractFloat,
             Pm::AbstractVector,
             Dm::AbstractVector,
             dm::AbstractVector,
             x0::AbstractFloat,
             y0::AbstractFloat,
             y1::AbstractFloat,
             ft::AbstractVector,
             FT::AbstractVector,
             mf::AbstractFloat,
             bPt::AbstractVector,
             bDt::AbstractVector,
             bdt::AbstractVector)
    return @. ((Dm-bDt)*mf*wd*y1^2+(((FT*bPt*ft-FT*ft*Pm)*wP*wd*x0+(bDt-Dm)*mf*wd)*y0+(bdt-dm)*mf^2*wd)*y1+((FT^2*ft^2*dm-FT^2*bdt*ft^2)*wP*wd*x0^2+(dm-bdt)*mf^2*wd)*y0+(FT^2*ft^2*Dm-FT^2*bDt*ft^2)*mf*wP*x0^2+(FT*bPt*ft-FT*ft*Pm)*mf^2*wP*x0)/((Dm-bDt)*mf*wd*y1^2+((FT*bPt*ft-FT*ft*Pm)*wP*wd*x0+(2*bDt-2*Dm)*mf*wd)*y0*y1+((FT*ft*Pm-FT*bPt*ft)*wP*wd*x0+(Dm-bDt)*mf*wd)*y0^2+(FT^2*ft^2*dm-FT^2*bdt*ft^2)*wP*wd*x0^2*y0+(FT^2*ft^2*Dm-FT^2*bDt*ft^2)*mf*wP*x0^2)
end
# for glass
function w2S(wd::AbstractFloat,
             Dm::AbstractVector,
             dm::AbstractVector,
             y0::AbstractFloat,
             mf::AbstractFloat,
             bDt::AbstractVector,
             bdt::AbstractVector)
    return @. ((dm-bdt)*wd*y0^2+((dm-bdt)*wd+(Dm-bDt)*mf)*y0+(Dm-bDt)*mf)/(wd*y0^2+mf^2)
end

# mineral
function SS(par::AbstractVector,
            wP::AbstractFloat,
            wd::AbstractFloat,
            bP::AbstractVector,
            bD::AbstractVector,
            bd::AbstractVector,
            dats::AbstractDict,
            channels::AbstractDict,
            anchors::AbstractDict,
            mf::Union{AbstractFloat,Nothing};
            ndrift::Integer=1,
            ndown::Integer=0,
            PAcutoff=nothing,
            verbose::Bool=false)
    drift = par[1:ndrift]
    down = vcat(0.0,par[ndrift+1:ndrift+ndown])
    mfrac = isnothing(mf) ? par[ndrift+ndown+1] : log(mf)
    adrift = isnothing(PAcutoff) ? drift : par[end-ndrift+1:end]
    out = 0.0
    for (refmat,dat) in dats
        (x0,y0,y1) = anchors[refmat]
        Pm,Dm,dm,ft,FT,mf,bPt,bDt,bdt =
            SSprep(bP,bD,bd,dat,channels,mfrac,drift,down;
                   PAcutoff=PAcutoff,adrift=adrift)
        out += SS(wP,wd,Pm,Dm,dm,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    end
    return out
end
function SS(wP::AbstractFloat,
            wd::AbstractFloat,
            Pm::AbstractVector,
            Dm::AbstractVector,
            dm::AbstractVector,
            x0::AbstractFloat,
            y0::AbstractFloat,
            y1::AbstractFloat,
            ft::AbstractVector,
            FT::AbstractVector,
            mf::AbstractFloat,
            bPt::AbstractVector,
            bDt::AbstractVector,
            bdt::AbstractVector)
    pred = predict(wP,wd,Pm,Dm,dm,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    ss = @. wP*(pred[:,"P"]-Pm)^2 + (pred[:,"D"]-Dm)^2 + wd*(pred[:,"d"]-dm)^2
    return sum(ss)
end
# glass
function SS(par::AbstractVector,
            wd::AbstractFloat,
            bD::AbstractVector,
            bd::AbstractVector,
            dats::AbstractDict,
            channels::AbstractDict,
            anchors::AbstractDict)
    mf = exp(par[1])
    out = 0.0
    for (refmat,dat) in dats
        y0 = anchors[refmat]
        Dm,dm,bDt,bdt = SSprep(bD,bd,dat,channels)
        out += SS(wd,Dm,dm,y0,mf,bDt,bdt)
    end
    return out
end
function SS(wd::AbstractFloat,
            Dm::AbstractVector,
            dm::AbstractVector,
            y0::AbstractFloat,
            mf::AbstractFloat,
            bDt::AbstractVector,
            bdt::AbstractVector)
    pred = predict(wd,Dm,dm,y0,mf,bDt,bdt)
    ss = @. (pred[:,"D"]-Dm)^2 + wd*(pred[:,"d"]-dm)^2
    return sum(ss)
end
export SS

# minerals
function SSprep(bP,bD,bd,dat,channels,mfrac,drift,down;
                PAcutoff=nothing,adrift=drift)
    t = dat.t
    T = dat.T 
    Pm = dat[:,channels["P"]]
    Dm = dat[:,channels["D"]]
    dm = dat[:,channels["d"]]
    ft = get_drift(Pm,t,drift;
                   PAcutoff=PAcutoff,adrift=adrift)
    FT = polyFac(down,T)
    mf = exp(mfrac)
    bPt = polyVal(bP,t)
    bDt = polyVal(bD,t)
    bdt = polyVal(bd,t)
    return Pm,Dm,dm,ft,FT,mf,bPt,bDt,bdt
end
# glass
function SSprep(bD,bd,dat,channels)
    t = dat.t
    Dm = dat[:,channels["D"]]
    dm = dat[:,channels["d"]]
    bDt = polyVal(bD,t)
    bdt = polyVal(bd,t)
    return Dm,dm,bDt,bdt
end

# minerals
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
        anchor = anchors[samp.group]
        return predict(dat,pars,blank,channels,anchor;debug=debug)
    end
end
function predict(dat::AbstractDataFrame,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 anchor::NamedTuple;
                 debug::Bool=false)
    bP = blank[:,channels["P"]]
    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    Pm,Dm,dm,ft,FT,mf,bPt,bDt,bdt =
        SSprep(bP,bD,bd,dat,channels,
               pars.mfrac,pars.drift,pars.down;
               PAcutoff=pars.PAcutoff,adrift=pars.adrift)
    return predict(pars.wP,pars.wd,Pm,Dm,dm,
                   anchor.x0,anchor.y0,anchor.y1,
                   ft,FT,mf,bPt,bDt,bdt)
end
function predict(wP::AbstractFloat,
                 wd::AbstractFloat,
                 Pm::AbstractVector,
                 Dm::AbstractVector,
                 dm::AbstractVector,
                 x0::AbstractFloat,
                 y0::AbstractFloat,
                 y1::AbstractFloat,
                 ft::AbstractVector,
                 FT::AbstractVector,
                 mf::AbstractFloat,
                 bPt::AbstractVector,
                 bDt::AbstractVector,
                 bdt::AbstractVector)
    S = w2S(wP,wd,Pm,Dm,dm,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    p = w2p(wP,wd,Pm,Dm,dm,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    return predict(S,p,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
end
function predict(S::AbstractVector,
                 p::AbstractVector,
                 x0::AbstractFloat,
                 y0::AbstractFloat,
                 y1::AbstractFloat,
                 ft::AbstractVector,
                 FT::AbstractVector,
                 mf::AbstractFloat,
                 bPt::AbstractVector,
                 bDt::AbstractVector,
                 bdt::AbstractVector)
    x = @. x0*(1-p)
    y = @. y1+(y0-y1)*p
    z = @. 1+x+y
    Pf = @. S*ft*FT*x/z + bPt
    Df = @. S*mf/z + bDt
    df = @. S*y/z + bdt
    return DataFrame(P=Pf,D=Df,d=df)
end
# glass
function predict(wd::AbstractFloat,
                 Dm::AbstractVector,
                 dm::AbstractVector,
                 y0::AbstractFloat,
                 mf::AbstractFloat,
                 bDt::AbstractVector,
                 bdt::AbstractVector)
    S = w2S(wd,Dm,dm,y0,mf,bDt,bdt)
    y = @. y0
    z = @. 1+y
    Df = @. S*mf/z + bDt
    df = @. S*y/z + bdt
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
export predict

function averat_jacobian(P,D,d,x,y)
    ns = length(P)
    z = @. 1 + x + y
    S = @. (D+x*P+y*d)*z/(1 + x^2 + y^2)
    dPdS = fill(x/z,ns)
    dDdS = fill(1/z,ns)
    dddS = fill(y/z,ns)
    dPdx = @. S/z - S*x/z^2
    dDdx = @. - S/z^2
    dddx = @. - S*y/z^2
    dPdy = @. - S*x/z^2
    dDdy = @. - S/z^2
    dddy = @. S/z - S*y/z^2
    J = zeros(ns+2,3*ns)
    J[1,1:ns] .= dPdx
    J[1,ns+1:2*ns] .= dDdx
    J[1,2*ns+1:3*ns] .= dddx
    J[2,1:ns] .= dPdy
    J[2,ns+1:2*ns] .= dDdy
    J[2,2*ns+1:3*ns] .= dddy
    J[3:end,1:ns] .= diagm(dPdS)
    J[3:end,ns+1:2*ns] .= diagm(dDdS)
    J[3:end,2*ns+1:3*ns] .= diagm(dddS)
    return J
end

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
