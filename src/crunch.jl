# for age standards
function getq(Pm::AbstractVector,
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
    return @. ((bdt-bDt+Pm+2*Dm)*mf^2*y1^2+((((-FT*bdt)+FT*bPt-2*FT*Pm-Dm*FT)*ft*mf^2*x0+((-bdt)+bDt-Pm-2*Dm)*mf^2)*y0+((-FT*dm)+FT*bdt+FT*bPt-2*FT*bDt-FT*Pm+2*Dm*FT)*ft*mf*x0+((-2*dm)+bdt-bDt-Pm)*mf)*y1+((2*FT^2*dm-FT^2*bdt+FT^2*bPt+Dm*FT^2)*ft^2*mf*x0^2+(2*FT*dm-2*FT*bdt+FT*bPt+FT*bDt-FT*Pm-Dm*FT)*ft*mf*x0+(2*dm-bdt+bDt+Pm)*mf)*y0+(FT^2*dm+FT^2*bPt-FT^2*bDt+2*Dm*FT^2)*ft^2*x0^2+((-FT*dm)+FT*bPt-FT*bDt-2*FT*Pm)*ft*x0)/((bdt-bDt+Pm+2*Dm)*mf^2*y1^2+((((-FT*bdt)+FT*bPt-2*FT*Pm-Dm*FT)*ft*mf^2*x0+((-2*bdt)+2*bDt-2*Pm-4*Dm)*mf^2)*y0+((-FT*dm)+FT*bdt+FT*bPt-2*FT*bDt-FT*Pm+2*Dm*FT)*ft*mf*x0)*y1+((FT*bdt-FT*bPt+2*FT*Pm+Dm*FT)*ft*mf^2*x0+(bdt-bDt+Pm+2*Dm)*mf^2)*y0^2+((2*FT^2*dm-FT^2*bdt+FT^2*bPt+Dm*FT^2)*ft^2*mf*x0^2+(FT*dm-FT*bdt-FT*bPt+2*FT*bDt+FT*Pm-2*Dm*FT)*ft*mf*x0)*y0+(FT^2*dm+FT^2*bPt-FT^2*bDt+2*Dm*FT^2)*ft^2*x0^2)
end

# mineral
function residuals(par::AbstractVector,
                   bP::AbstractVector,
                   bD::AbstractVector,
                   bd::AbstractVector,
                   dats::AbstractDict,
                   channels::AbstractDict,
                   anchors::AbstractDict,
                   mf::Union{AbstractFloat,Nothing};
                   ndrift::Integer=1,
                   ndown::Integer=0,
                   PAcutoff=nothing)
    drift = par[1:ndrift]
    down = vcat(0.0,par[ndrift+1:ndrift+ndown])
    mfrac = isnothing(mf) ? par[ndrift+ndown+1] : log(mf)
    adrift = isnothing(PAcutoff) ? drift : par[end-ndrift+1:end]
    dP = Float64[]
    dD = Float64[]
    dd = Float64[]
    for (refmat,dat) in dats
        (x0,y0,y1) = anchors[refmat]
        Pm,Dm,dm,ft,FT,mf,bPt,bDt,bdt =
            SSprep(bP,bD,bd,dat,channels,mfrac,drift,down;
                   PAcutoff=PAcutoff,adrift=adrift)
        dPnew, dDnew, ddnew = residuals(Pm,Dm,dm,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
        append!(dP,dPnew)
        append!(dD,dDnew)
        append!(dd,ddnew)
    end
    return dP,dD,dd
end
function residuals(Pm::AbstractVector,
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
    pred = predict(Pm,Dm,dm,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    dP = pred[:,"P"] .- Pm
    dD = pred[:,"D"] .- Dm
    dd = pred[:,"d"] .- dm
    return dP, dD, dd
end
# glass
function residuals(par::AbstractVector,
                   bD::AbstractVector,
                   bd::AbstractVector,
                   dats::AbstractDict,
                   channels::AbstractDict,
                   anchors::AbstractDict)
    mf = exp(par[1])
    out = 0.0
    dD = Float64[]
    dd = Float64[]
    for (refmat,dat) in dats
        y0 = anchors[refmat]
        Dm,dm,bDt,bdt = SSprep(bD,bd,dat,channels)
        dDnew, ddnew = residuals(Dm,dm,y0,mf,bDt,bdt)
        append!(dD,dDnew)
        append!(dd,ddnew)
    end
    return dD, dd
end
function residuals(Dm::AbstractVector,
                   dm::AbstractVector,
                   y0::AbstractFloat,
                   mf::AbstractFloat,
                   bDt::AbstractVector,
                   bdt::AbstractVector)
    pred = predict(Dm,dm,y0,mf,bDt,bdt)
    dD = pred[:,"D"] .- Dm
    dd = pred[:,"d"] .- dm
    return dD, dd
end

# mineral
function SS(par::AbstractVector,
            bP::AbstractVector,
            bD::AbstractVector,
            bd::AbstractVector,
            dats::AbstractDict,
            channels::AbstractDict,
            anchors::AbstractDict,
            mf::Union{AbstractFloat,Nothing};
            ndrift::Integer=1,
            ndown::Integer=0,
            PAcutoff=nothing)
    dP, dD, dd = residuals(par,bP,bD,bd,dats,channels,anchors,mf;
                           ndrift=ndrift,ndown=ndown,PAcutoff=PAcutoff)
    return sum(dP.^2) + sum(dD.^2) + sum(dd.^2)
end
function SS(Pm::AbstractVector,
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
    dP, dD, dd = residuals(Pm,Dm,dm,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    return sum(dP.^2) + sum(dD.^2) + sum(dd.^2)
end
# glass
function SS(par::AbstractVector,
            bD::AbstractVector,
            bd::AbstractVector,
            dats::AbstractDict,
            channels::AbstractDict,
            anchors::AbstractDict)
    dD, dd = residuals(par,bD,bd,dats,channels,anchors)
    return sum(dD.^2) + sum(dd.^2)
end
function SS(Dm::AbstractVector,
            dm::AbstractVector,
            y0::AbstractFloat,
            mf::AbstractFloat,
            bDt::AbstractVector,
            bdt::AbstractVector)
    dD, dd = residuals(Dm,dm,y0,mf,bDt,bdt)
    return sum(dD.^2) + sum(dd.^2)
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
    return predict(Pm,Dm,dm,
                   anchor.x0,anchor.y0,anchor.y1,
                   ft,FT,mf,bPt,bDt,bdt)
end
function predict(Pm::AbstractVector,
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
    S = @. Pm + Dm + dm
    q = getq(Pm,Dm,dm,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    return predict(S,q,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
end
function predict(S::AbstractVector,
                 q::AbstractVector,
                 x0::AbstractFloat,
                 y0::AbstractFloat,
                 y1::AbstractFloat,
                 ft::AbstractVector,
                 FT::AbstractVector,
                 mf::AbstractFloat,
                 bPt::AbstractVector,
                 bDt::AbstractVector,
                 bdt::AbstractVector)
    x = @. x0*(1-q)*ft*FT
    y = @. (y1+(y0-y1)*q)*mf
    z = @. 1+x+y
    Pf = @. S*x/z + bPt
    Df = @. S/z + bDt
    df = @. S*y/z + bdt
    return DataFrame(P=Pf,D=Df,d=df)
end
# glass
function predict(Dm::AbstractVector,
                 dm::AbstractVector,
                 y0::AbstractFloat,
                 mf::AbstractFloat,
                 bDt::AbstractVector,
                 bdt::AbstractVector)
    S = @. Dm + dm
    y = @. y0*mf
    z = @. 1+y
    Df = @. S/z + bDt
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
