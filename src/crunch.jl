# for age standards
function getS(Pm,Dm,dm,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    S = @. (((dm-bdt)*mf^2+(Dm-bDt)*mf)*y1^2+((((bDt-Dm)*mf+(FT*bPt-FT*Pm)*ft)*x0+(2*bdt-2*dm)*mf^2+(2*bDt-2*Dm)*mf)*y0+((FT*Pm-FT*bPt)*ft+dm-bdt)*mf^2*x0)*y1+(((FT^2*dm-FT^2*bdt)*ft^2+(FT*Pm-FT*bPt)*ft)*x0^2+((Dm-bDt)*mf+(FT*Pm-FT*bPt)*ft)*x0+(dm-bdt)*mf^2+(Dm-bDt)*mf)*y0^2+(((Dm*FT^2-FT^2*bDt)*ft^2*mf+(FT^2*dm-FT^2*bdt)*ft^2)*x0^2+((FT*bPt-FT*Pm)*ft-dm+bdt)*mf^2*x0)*y0+((FT*Pm-FT*bPt)*ft*mf^2+(Dm*FT^2-FT^2*bDt)*ft^2*mf)*x0^2)/(mf^2*y1^2-2*mf^2*y0*y1+(FT^2*ft^2*x0^2+mf^2)*y0^2+FT^2*ft^2*mf^2*x0^2)
    return S
end
# for glass
function getS(Dm,dm,y0,mf,bDt,bdt)
    S = @. ((dm-bdt)*y0^2+((Dm-bDt)*mf+dm-bdt)*y0+(Dm-bDt)*mf)/(y0^2+mf^2)
    return S
end
export getS
function getp(Pm,Dm,dm,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    p = @. ((bDt-Dm)*mf*y1^2+(((FT*Pm-FT*bPt)*ft*x0+(Dm-bDt)*mf)*y0+(dm-bdt)*mf^2)*y1+((FT^2*bdt-FT^2*dm)*ft^2*x0^2+(bdt-dm)*mf^2)*y0+(FT^2*bDt-Dm*FT^2)*ft^2*mf*x0^2+(FT*Pm-FT*bPt)*ft*mf^2*x0)/((bDt-Dm)*mf*y1^2+((FT*Pm-FT*bPt)*ft*x0+(2*Dm-2*bDt)*mf)*y0*y1+((FT*bPt-FT*Pm)*ft*x0+(bDt-Dm)*mf)*y0^2+(FT^2*bdt-FT^2*dm)*ft^2*x0^2*y0+(FT^2*bDt-Dm*FT^2)*ft^2*mf*x0^2)
    return p
end
export getp

# isotopic ratios in matrix matched mineral standards
function SS(t,T,Pm,Dm,dm,x0,y0,y1,drift,down,mfrac,bP,bD,bd;
            PAcutoff=nothing,adrift=drift)
    pred = predict(t,T,Pm,Dm,dm,x0,y0,y1,drift,down,mfrac,bP,bD,bd;
                   PAcutoff=PAcutoff,adrift=adrift)
    S = @. (pred[:,"P"]-Pm)^2 + (pred[:,"D"]-Dm)^2 + (pred[:,"d"]-dm)^2
    return sum(S)
end
# isotopic ratios in glass
function SS(t,Dm,dm,y0,mfrac,bD,bd)
    pred = predict(t,Dm,dm,y0,mfrac,bD,bd)
    S = @. (pred[:,"D"]-Dm)^2 + (pred[:,"d"]-dm)^2
    return sum(S)
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
function get_drift(Pm::AbstractVector,
                   t::AbstractVector,
                   pars::NamedTuple)
    return get_drift(Pm,t,pars.drift;
                     PAcutoff=pars.PAcutoff,
                     adrift=pars.adrift)
end

# isotopic ratios
function predict(t,T,Pm,Dm,dm,x0,y0,y1,drift,down,mfrac,bP,bD,bd;
                 PAcutoff=nothing,adrift=drift,debug::Bool=false)
    ft = get_drift(Pm,t,drift;PAcutoff=PAcutoff,adrift=adrift)
    FT = polyFac(down,T)
    mf = exp(mfrac)
    bPt = polyVal(bP,t)
    bDt = polyVal(bD,t)
    bdt = polyVal(bd,t)
    S = getS(Pm,Dm,dm,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    p = getp(Pm,Dm,dm,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    x = @. x0*(1-p)
    y = @. y1+(y0-y1)*p
    z = @. 1+x+y
    Pf = @. S*ft*FT*x/z + bPt
    Df = @. S*mf/z + bDt
    df = @. S*y/z + bdt
    if debug
        @infiltrate
    end
    return DataFrame(P=Pf,D=Df,d=df)
end
# isotopic ratios for glass
function predict(t,Dm,dm,y0,mfrac,bD,bd;debug::Bool=false)
    mf = exp(mfrac)
    bDt = polyVal(bD,t)
    bdt = polyVal(bd,t)
    S = getS(Dm,dm,y0,mf,bDt,bdt)
    x = 0
    y = y0
    z = 1+x+y
    Df = @. S*mf/z + bDt
    df = @. S*y/z + bdt
    return DataFrame(D=Df,d=df)
end
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
# minerals
function predict(dat::AbstractDataFrame,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 anchor::NamedTuple;
                 debug::Bool=false)
    t = dat.t
    T = dat.T
    Pm = dat[:,channels["P"]]
    Dm = dat[:,channels["D"]]
    dm = dat[:,channels["d"]]
    bP = blank[:,channels["P"]]
    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    return predict(t,T,Pm,Dm,dm,
                   anchor.x0,anchor.y0,anchor.y1,
                   pars.drift,pars.down,pars.mfrac,bP,bD,bd;
                   PAcutoff=pars.PAcutoff,adrift=pars.adrift,
                   debug=debug)
end
# glass
function predict(dat::AbstractDataFrame,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 y0::AbstractFloat;
                 debug::Bool=false)
    t = dat.t
    Dm = dat[:,channels["D"]]
    dm = dat[:,channels["d"]]
    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    return predict(t,Dm,dm,y0,pars.mfrac,bD,bd;debug=debug)
end
# concentrations
function predict(samp::Sample,
                 ef::AbstractVector,
                 blank::AbstractDataFrame,
                 elements::AbstractDataFrame,
                 internal::AbstractString;
                 debug::Bool=debug)
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
