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
    return @. ((FT*(bDt-Dm)*ft*mf*vP*x0*y0+(FT*dm-FT*bdt)*ft*vP*x0)*y1+(mf^2*(Pm*vD-bPt*vD)*x0^2+FT*(Dm-bDt)*ft*mf*vP*x0)*y0^2+(FT*bdt-FT*dm)*ft*vP*x0*y0+(Pm-bPt)*vd*x0^2)/(FT^2*ft^2*vP*y1^2-2*FT^2*ft^2*vP*y0*y1+(mf^2*vD*x0^2+FT^2*ft^2*vP)*y0^2+vd*x0^2)
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
    return @. -(FT^2*(bDt-Dm)*ft^2*vP*y1^2+(ft*mf*(FT*Pm*vD-FT*bPt*vD)*x0+FT^2*(2*Dm-2*bDt)*ft^2*vP)*y0*y1+(ft*mf*(FT*bPt*vD-FT*Pm*vD)*x0+FT^2*(bDt-Dm)*ft^2*vP)*y0^2+mf*(bdt*vD-dm*vD)*x0^2*y0+(bDt-Dm)*vd*x0^2)/(FT^2*ft^2*vP*y1^2-2*FT^2*ft^2*vP*y0*y1+(mf^2*vD*x0^2+FT^2*ft^2*vP)*y0^2+vd*x0^2)
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
    return @. ((dm-bdt)*mf*vD*y0+(Dm-bDt)*vd)/(mf^2*vD*y0^2+vd)
end

# mineral
function LL(par::AbstractVector,
            bP::AbstractVector,
            bD::AbstractVector,
            bd::AbstractVector,
            dats::AbstractDict,
            dt::AbstractDict,
            channels::AbstractDict,
            anchors::AbstractDict,
            mf::Union{AbstractFloat,Nothing};
            ndrift::Integer=1,
            ndown::Integer=0,
            PAcutoff::Union{AbstractFloat,Nothing}=nothing,
            dead::AbstractFloat=0.0)
    drift = par[1:ndrift]
    down = vcat(0.0,par[ndrift+1:ndrift+ndown])
    mfrac = isnothing(mf) ? par[ndrift+ndown+1] : log(mf)
    adrift = isnothing(PAcutoff) ? drift : par[end-ndrift+1:end]
    out = 0.0
    for (refmat,dat) in dats
        (x0,y0,y1) = anchors[refmat]
        Pm,Dm,dm,vP,vD,vd,ft,FT,mf,bPt,bDt,bdt =
            LLprep(bP,bD,bd,dat,dt,channels,mfrac,drift,down;
                   PAcutoff=PAcutoff,adrift=adrift,dead=dead)
        out += LL(Pm,Dm,dm,vP,vD,vd,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    end
    return out
end
function LL(Pm::AbstractVector,
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
    return sum(@. log(vP*vD*vd) + (dP^2)/vP + (dD^2)/vD + (dd^2)/vd )/2 
end
# glass
function LL(par::AbstractVector,
            bD::AbstractVector,
            bd::AbstractVector,
            dats::AbstractDict,
            dt::AbstractDict,
            channels::AbstractDict,
            anchors::AbstractDict;
            dead::AbstractFloat=0.0)
    mf = exp(par[1])
    out = 0.0
    for (refmat,dat) in dats
        y0 = anchors[refmat]
        Dm,dm,vD,vd,bDt,bdt = LLprep(bD,bd,dat,dt,channels;dead=dead)
        out += LL(Dm,dm,vD,vd,y0,mf,bDt,bdt)
    end
    return out
end
function LL(Dm::AbstractVector,
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
    return sum(@. log(vD*vd) + (dD^2)/vD + (dd^2)/vd )/2
end
export LL

# minerals
function LLprep(bP::AbstractVector,
                bD::AbstractVector,
                bd::AbstractVector,
                dat::AbstractDataFrame,
                dt::AbstractDict,
                channels::AbstractDict,
                mfrac::AbstractFloat,
                drift::AbstractVector,
                down::AbstractVector;
                PAcutoff::Union{AbstractFloat,Nothing}=nothing,
                adrift::AbstractVector=drift,
                dead::AbstractFloat=0.0)
    t = dat.t
    T = dat.T 
    Pm = dat[:,channels["P"]]
    Dm = dat[:,channels["D"]]
    dm = dat[:,channels["d"]]
    vP = var_cps(Pm,dt[channels["P"]],dead)
    vD = var_cps(Dm,dt[channels["D"]],dead)
    vd = var_cps(dm,dt[channels["d"]],dead)
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
function LLprep(bD::AbstractVector,
                bd::AbstractVector,
                dat::AbstractDataFrame,
                dt::AbstractDict,
                channels::AbstractDict;
                dead::AbstractFloat=0.0)
    t = dat.t
    Dm = dat[:,channels["D"]]
    dm = dat[:,channels["d"]]
    vD = var_cps(Dm,dt[channels["D"]],dead)
    vd = var_cps(dm,dt[channels["d"]],dead)
    bDt = polyVal(bD,t)
    bdt = polyVal(bd,t)
    return Dm,dm,vD,vd,bDt,bdt
end

# minerals
function predict(samp::Sample,
                 dt::AbstractDict,
                 method::AbstractString,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 standards::AbstractDict,
                 glass::AbstractDict;
                 debug::Bool=false)
    anchors = getAnchors(method,standards,glass)
    return predict(samp,dt,pars,blank,channels,anchors;debug=debug)
end
function predict(samp::Sample,
                 dt::AbstractDict,
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
        return predict(dat,dt,pars,blank,channels,anchor;debug=debug)
    end
end
function predict(dat::AbstractDataFrame,
                 dt::AbstractDict,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 anchor::NamedTuple;
                 debug::Bool=false)
    bP = blank[:,channels["P"]]
    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    Pm,Dm,dm,vP,vD,vd,ft,FT,mf,bPt,bDt,bdt =
        LLprep(bP,bD,bd,dat,dt,channels,
               pars.mfrac,pars.drift,pars.down;
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
    df = @. D*mf*y0 + P*ft*FT*(y1-y0)/x0 + bdt
    return DataFrame(P=Pf,D=Df,d=df)
end
# glass
function predict(dat::AbstractDataFrame,
                 dt::AbstractDict,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 y0::AbstractFloat;
                 debug::Bool=false)
    mf = exp(pars.mfrac)
    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    Dm,dm,vD,vd,bDt,bdt = LLprep(bD,bd,dat,dt,channels)
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
