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
                 PAcutoff=nothing,adrift=drift)
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
    return DataFrame(P=Pf,D=Df,d=df)
end
# isotopic ratios for glass
function predict(t,Dm,dm,y0,mfrac,bD,bd)
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
                 glass::AbstractDict)
    anchors = getAnchors(method,standards,glass)
    return predict(samp,pars,blank,channels,anchors)
end
function predict(samp::Sample,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 anchors::AbstractDict)
    if samp.group == "sample"
        PTerror("notStandard")
    else
        dat = windowData(samp;signal=true)
        anchor = anchors[samp.group]
        return predict(dat,pars,blank,channels,anchor)
    end
end
# minerals
function predict(dat::AbstractDataFrame,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 anchor::NamedTuple)
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
                   PAcutoff=pars.PAcutoff,adrift=pars.adrift)
end
# glass
function predict(dat::AbstractDataFrame,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 y0::AbstractFloat)
    t = dat.t
    Dm = dat[:,channels["D"]]
    dm = dat[:,channels["d"]]
    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    return predict(t,Dm,dm,y0,pars.mfrac,bD,bd)
end
# concentrations
function predict(samp::Sample,
                 ef::AbstractVector,
                 blank::AbstractDataFrame,
                 elements::AbstractDataFrame,
                 internal::AbstractString)
    if samp.group in collect(keys(_PT["glass"]))
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
        PTerror("notStandard")
    end
end
export predict

function get_covmat_averat(P,D,d,x,y)
    S = @. (d*y^2+((d+P)*x+d+D)*y+P*x^2+(P+D)*x+D)/(y^2+x^2+1)
    z = @. 1 + x + y
    dSSdS = @. ((2*y*(S*y/z-d)) + (2*x*(S*x/z-P)) + (2*(D-S/z)))/z
    dSSdx = sum( @. (2*S*y*(d-(S*y)/z))/z^2 + 2*((S*x)/z^2-S/z)*(P-(S*x)/z) + (2*S*(D-S/z))/z^2 )
    dSSdy = sum( @. (2*S*x*(P-(S*x)/z))/z^2 + 2*((S*y)/z^2-S/z)*(d-(S*y)/z) + (2*S*(D-S/z))/z^2 )
    d2SSdS2 = @. 2*(y^2+x^2+1)/z^2
    d2SSdSdx = @. (2*y*(d-(S*y)/z))/z^2-(2*(P-(S*x)/z))/z+(2*x*(P-(S*x)/z))/z^2-(2*x*((S*x)/z^2-S/z))/z+(2*(D-S/z))/z^2-(2*S*y^2)/z^3-(2*S)/z^3
    d2SSdSdy = @. (-(2*(d-(S*y)/z))/z)+(2*y*(d-(S*y)/z))/z^2+(2*x*(P-(S*x)/z))/z^2-(2*y*((S*y)/z^2-S/z))/z+(2*(D-S/z))/z^2-(2*S*x^2)/z^3-(2*S)/z^3
    d2SSdx2 = sum( @. 2*((S*x)/z^2-S/z)^2-(4*S*y*(d-(S*y)/z))/z^3+2*((2*S)/z^2-(2*S*x)/z^3)*(P-(S*x)/z)-(4*S*(D-S/z))/z^3+(2*S^2*y^2)/z^4+(2*S^2)/z^4 )
    d2SSdy2 = sum( @. 2*((S*y)/z^2-S/z)^2+2*((2*S)/z^2-(2*S*y)/z^3)*(d-(S*y)/z)-(4*S*x*(P-(S*x)/z))/z^3-(4*S*(D-S/z))/z^3+(2*S^2*x^2)/z^4+(2*S^2)/z^4 )
    d2SSdxdy = sum( @. (2*S*(d-(S*y)/z))/z^2-(4*S*y*(d-(S*y)/z))/z^3+2*(S/z^2-(2*S*x)/z^3)*(P-(S*x)/z)+(2*S*y*((S*y)/z^2-S/z))/z^2+(2*S*x*((S*x)/z^2-S/z))/z^2-(4*S*(D-S/z))/z^3+(2*S^2)/z^4 )
    ns = length(P)
    J = zeros(1,ns+2)
    J[1,1:ns] .= dSSdS
    J[1,ns+1] = dSSdx
    J[1,ns+2] = dSSdy
    H = zeros(ns+2,ns+2)
    H[diagind(H)[1:ns]] .= d2SSdS2
    H[ns+1,1:ns] .= d2SSdSdx
    H[ns+2,1:ns] .= d2SSdSdy
    H[ns+1,ns+1] = d2SSdx2
    H[ns+2,ns+2] = d2SSdy2
    H[ns+1,ns+2] = d2SSdxdy
    H[ns+2,ns+1] = d2SSdxdy
    covmat = inv(-transpose(J)*J*H)
    return covmat
end
