# minerals
function getP(Pm::AbstractVector,
              Dm::AbstractVector,
              dm::AbstractVector,
              O::Matrix,
              x0::AbstractFloat,
              y0::AbstractFloat,
              y1::AbstractFloat,
              ft::AbstractVector,
              FT::AbstractVector,
              mf::AbstractFloat,
              bPt::AbstractVector,
              bDt::AbstractVector,
              bdt::AbstractVector)
    return @. -(((((2*FT*O[2,1]+2*FT*O[1,2])*O[3,3]+((-FT*O[3,1])-FT*O[1,3])*O[3,2]-FT*O[2,3]*O[3,1]-FT*O[1,3]*O[2,3])*Pm+4*Dm*FT*O[2,2]*O[3,3]-Dm*FT*O[3,2]^2-2*Dm*FT*O[2,3]*O[3,2]-Dm*FT*O[2,3]^2)*ft*mf*x0*y0+(((-4*FT*O[2,2]*O[3,3])+FT*O[3,2]^2+2*FT*O[2,3]*O[3,2]+FT*O[2,3]^2)*dm+((FT*O[2,1]+FT*O[1,2])*O[3,2]-2*FT*O[2,2]*O[3,1]+(FT*O[2,1]+FT*O[1,2])*O[2,3]-2*FT*O[1,3]*O[2,2])*Pm)*ft*x0)*y1+((((-4*O[1,1]*O[3,3])+O[3,1]^2+2*O[1,3]*O[3,1]+O[1,3]^2)*Pm+((-2*Dm*O[2,1])-2*Dm*O[1,2])*O[3,3]+(Dm*O[3,1]+Dm*O[1,3])*O[3,2]+Dm*O[2,3]*O[3,1]+Dm*O[1,3]*O[2,3])*mf^2*x0^2+((((-2*FT*O[2,1])-2*FT*O[1,2])*O[3,3]+(FT*O[3,1]+FT*O[1,3])*O[3,2]+FT*O[2,3]*O[3,1]+FT*O[1,3]*O[2,3])*Pm-4*Dm*FT*O[2,2]*O[3,3]+Dm*FT*O[3,2]^2+2*Dm*FT*O[2,3]*O[3,2]+Dm*FT*O[2,3]^2)*ft*mf*x0)*y0^2+((((2*O[2,1]+2*O[1,2])*O[3,3]+((-O[3,1])-O[1,3])*O[3,2]-O[2,3]*O[3,1]-O[1,3]*O[2,3])*dm+((-4*O[1,1]*O[3,2])+(2*O[2,1]+2*O[1,2])*O[3,1]-4*O[1,1]*O[2,3]+2*O[1,3]*O[2,1]+2*O[1,2]*O[1,3])*Pm+((-Dm*O[2,1])-Dm*O[1,2])*O[3,2]+2*Dm*O[2,2]*O[3,1]+((-Dm*O[2,1])-Dm*O[1,2])*O[2,3]+2*Dm*O[1,3]*O[2,2])*mf*x0^2+((4*FT*O[2,2]*O[3,3]-FT*O[3,2]^2-2*FT*O[2,3]*O[3,2]-FT*O[2,3]^2)*dm+(((-FT*O[2,1])-FT*O[1,2])*O[3,2]+2*FT*O[2,2]*O[3,1]+((-FT*O[2,1])-FT*O[1,2])*O[2,3]+2*FT*O[1,3]*O[2,2])*Pm)*ft*x0)*y0+(((O[2,1]+O[1,2])*O[3,2]-2*O[2,2]*O[3,1]+(O[2,1]+O[1,2])*O[2,3]-2*O[1,3]*O[2,2])*dm+((-4*O[1,1]*O[2,2])+O[2,1]^2+2*O[1,2]*O[2,1]+O[1,2]^2)*Pm)*x0^2)/((4*FT^2*O[2,2]*O[3,3]-FT^2*O[3,2]^2-2*FT^2*O[2,3]*O[3,2]-FT^2*O[2,3]^2)*ft^2*y1^2+(((((-4*FT*O[2,1])-4*FT*O[1,2])*O[3,3]+(2*FT*O[3,1]+2*FT*O[1,3])*O[3,2]+2*FT*O[2,3]*O[3,1]+2*FT*O[1,3]*O[2,3])*ft*mf*x0+((-8*FT^2*O[2,2]*O[3,3])+2*FT^2*O[3,2]^2+4*FT^2*O[2,3]*O[3,2]+2*FT^2*O[2,3]^2)*ft^2)*y0+(((-2*FT*O[2,1])-2*FT*O[1,2])*O[3,2]+4*FT*O[2,2]*O[3,1]+((-2*FT*O[2,1])-2*FT*O[1,2])*O[2,3]+4*FT*O[1,3]*O[2,2])*ft*x0)*y1+((4*O[1,1]*O[3,3]-O[3,1]^2-2*O[1,3]*O[3,1]-O[1,3]^2)*mf^2*x0^2+((4*FT*O[2,1]+4*FT*O[1,2])*O[3,3]+((-2*FT*O[3,1])-2*FT*O[1,3])*O[3,2]-2*FT*O[2,3]*O[3,1]-2*FT*O[1,3]*O[2,3])*ft*mf*x0+(4*FT^2*O[2,2]*O[3,3]-FT^2*O[3,2]^2-2*FT^2*O[2,3]*O[3,2]-FT^2*O[2,3]^2)*ft^2)*y0^2+((4*O[1,1]*O[3,2]+((-2*O[2,1])-2*O[1,2])*O[3,1]+4*O[1,1]*O[2,3]-2*O[1,3]*O[2,1]-2*O[1,2]*O[1,3])*mf*x0^2+((2*FT*O[2,1]+2*FT*O[1,2])*O[3,2]-4*FT*O[2,2]*O[3,1]+(2*FT*O[2,1]+2*FT*O[1,2])*O[2,3]-4*FT*O[1,3]*O[2,2])*ft*x0)*y0+(4*O[1,1]*O[2,2]-O[2,1]^2-2*O[1,2]*O[2,1]-O[1,2]^2)*x0^2)
end
function getD(Pm::AbstractVector,
              Dm::AbstractVector,
              dm::AbstractVector,
              O::Matrix,
              x0::AbstractFloat,
              y0::AbstractFloat,
              y1::AbstractFloat,
              ft::AbstractVector,
              FT::AbstractVector,
              mf::AbstractFloat,
              bPt::AbstractVector,
              bDt::AbstractVector,
              bdt::AbstractVector)
    return @. ((((2*FT^2*O[2,1]+2*FT^2*O[1,2])*O[3,3]+((-FT^2*O[3,1])-FT^2*O[1,3])*O[3,2]-FT^2*O[2,3]*O[3,1]-FT^2*O[1,3]*O[2,3])*Pm+4*Dm*FT^2*O[2,2]*O[3,3]-Dm*FT^2*O[3,2]^2-2*Dm*FT^2*O[2,3]*O[3,2]-Dm*FT^2*O[2,3]^2)*ft^2*y1^2+(((((-4*FT*O[1,1]*O[3,3])+FT*O[3,1]^2+2*FT*O[1,3]*O[3,1]+FT*O[1,3]^2)*Pm+((-2*Dm*FT*O[2,1])-2*Dm*FT*O[1,2])*O[3,3]+(Dm*FT*O[3,1]+Dm*FT*O[1,3])*O[3,2]+Dm*FT*O[2,3]*O[3,1]+Dm*FT*O[1,3]*O[2,3])*ft*mf*x0+((((-4*FT^2*O[2,1])-4*FT^2*O[1,2])*O[3,3]+(2*FT^2*O[3,1]+2*FT^2*O[1,3])*O[3,2]+2*FT^2*O[2,3]*O[3,1]+2*FT^2*O[1,3]*O[2,3])*Pm-8*Dm*FT^2*O[2,2]*O[3,3]+2*Dm*FT^2*O[3,2]^2+4*Dm*FT^2*O[2,3]*O[3,2]+2*Dm*FT^2*O[2,3]^2)*ft^2)*y0+((((-2*FT*O[2,1])-2*FT*O[1,2])*O[3,3]+(FT*O[3,1]+FT*O[1,3])*O[3,2]+FT*O[2,3]*O[3,1]+FT*O[1,3]*O[2,3])*dm+((-2*FT*O[1,1]*O[3,2])+(FT*O[2,1]+FT*O[1,2])*O[3,1]-2*FT*O[1,1]*O[2,3]+FT*O[1,3]*O[2,1]+FT*O[1,2]*O[1,3])*Pm+((-2*Dm*FT*O[2,1])-2*Dm*FT*O[1,2])*O[3,2]+4*Dm*FT*O[2,2]*O[3,1]+((-2*Dm*FT*O[2,1])-2*Dm*FT*O[1,2])*O[2,3]+4*Dm*FT*O[1,3]*O[2,2])*ft*x0)*y1+(((4*FT*O[1,1]*O[3,3]-FT*O[3,1]^2-2*FT*O[1,3]*O[3,1]-FT*O[1,3]^2)*Pm+(2*Dm*FT*O[2,1]+2*Dm*FT*O[1,2])*O[3,3]+((-Dm*FT*O[3,1])-Dm*FT*O[1,3])*O[3,2]-Dm*FT*O[2,3]*O[3,1]-Dm*FT*O[1,3]*O[2,3])*ft*mf*x0+(((2*FT^2*O[2,1]+2*FT^2*O[1,2])*O[3,3]+((-FT^2*O[3,1])-FT^2*O[1,3])*O[3,2]-FT^2*O[2,3]*O[3,1]-FT^2*O[1,3]*O[2,3])*Pm+4*Dm*FT^2*O[2,2]*O[3,3]-Dm*FT^2*O[3,2]^2-2*Dm*FT^2*O[2,3]*O[3,2]-Dm*FT^2*O[2,3]^2)*ft^2)*y0^2+(((4*O[1,1]*O[3,3]-O[3,1]^2-2*O[1,3]*O[3,1]-O[1,3]^2)*dm+2*Dm*O[1,1]*O[3,2]+((-Dm*O[2,1])-Dm*O[1,2])*O[3,1]+2*Dm*O[1,1]*O[2,3]-Dm*O[1,3]*O[2,1]-Dm*O[1,2]*O[1,3])*mf*x0^2+(((2*FT*O[2,1]+2*FT*O[1,2])*O[3,3]+((-FT*O[3,1])-FT*O[1,3])*O[3,2]-FT*O[2,3]*O[3,1]-FT*O[1,3]*O[2,3])*dm+(2*FT*O[1,1]*O[3,2]+((-FT*O[2,1])-FT*O[1,2])*O[3,1]+2*FT*O[1,1]*O[2,3]-FT*O[1,3]*O[2,1]-FT*O[1,2]*O[1,3])*Pm+(2*Dm*FT*O[2,1]+2*Dm*FT*O[1,2])*O[3,2]-4*Dm*FT*O[2,2]*O[3,1]+(2*Dm*FT*O[2,1]+2*Dm*FT*O[1,2])*O[2,3]-4*Dm*FT*O[1,3]*O[2,2])*ft*x0)*y0+((2*O[1,1]*O[3,2]+((-O[2,1])-O[1,2])*O[3,1]+2*O[1,1]*O[2,3]-O[1,3]*O[2,1]-O[1,2]*O[1,3])*dm+4*Dm*O[1,1]*O[2,2]-Dm*O[2,1]^2-2*Dm*O[1,2]*O[2,1]-Dm*O[1,2]^2)*x0^2)/((4*FT^2*O[2,2]*O[3,3]-FT^2*O[3,2]^2-2*FT^2*O[2,3]*O[3,2]-FT^2*O[2,3]^2)*ft^2*y1^2+(((((-4*FT*O[2,1])-4*FT*O[1,2])*O[3,3]+(2*FT*O[3,1]+2*FT*O[1,3])*O[3,2]+2*FT*O[2,3]*O[3,1]+2*FT*O[1,3]*O[2,3])*ft*mf*x0+((-8*FT^2*O[2,2]*O[3,3])+2*FT^2*O[3,2]^2+4*FT^2*O[2,3]*O[3,2]+2*FT^2*O[2,3]^2)*ft^2)*y0+(((-2*FT*O[2,1])-2*FT*O[1,2])*O[3,2]+4*FT*O[2,2]*O[3,1]+((-2*FT*O[2,1])-2*FT*O[1,2])*O[2,3]+4*FT*O[1,3]*O[2,2])*ft*x0)*y1+((4*O[1,1]*O[3,3]-O[3,1]^2-2*O[1,3]*O[3,1]-O[1,3]^2)*mf^2*x0^2+((4*FT*O[2,1]+4*FT*O[1,2])*O[3,3]+((-2*FT*O[3,1])-2*FT*O[1,3])*O[3,2]-2*FT*O[2,3]*O[3,1]-2*FT*O[1,3]*O[2,3])*ft*mf*x0+(4*FT^2*O[2,2]*O[3,3]-FT^2*O[3,2]^2-2*FT^2*O[2,3]*O[3,2]-FT^2*O[2,3]^2)*ft^2)*y0^2+((4*O[1,1]*O[3,2]+((-2*O[2,1])-2*O[1,2])*O[3,1]+4*O[1,1]*O[2,3]-2*O[1,3]*O[2,1]-2*O[1,2]*O[1,3])*mf*x0^2+((2*FT*O[2,1]+2*FT*O[1,2])*O[3,2]-4*FT*O[2,2]*O[3,1]+(2*FT*O[2,1]+2*FT*O[1,2])*O[2,3]-4*FT*O[1,3]*O[2,2])*ft*x0)*y0+(4*O[1,1]*O[2,2]-O[2,1]^2-2*O[1,2]*O[2,1]-O[1,2]^2)*x0^2)
end
# glass
function getD(Dm::AbstractVector,
              dm::AbstractVector,
              O::Matrix,
              y0::AbstractFloat,
              mf::AbstractFloat,
              bDt::AbstractVector,
              bdt::AbstractVector)
    return @. ((2*O[2,2]*dm+Dm*O[2,1]+Dm*O[1,2])*mf*y0+(O[2,1]+O[1,2])*dm+2*Dm*O[1,1])/(2*O[2,2]*mf^2*y0^2+(2*O[2,1]+2*O[1,2])*mf*y0+2*O[1,1])
end

# mineral
function LL(par::AbstractVector,
            bP::AbstractVector,
            bD::AbstractVector,
            bd::AbstractVector,
            dats::AbstractDict,
            channels::AbstractDict,
            anchors::AbstractDict,
            mf::Union{AbstractFloat,Nothing};
            ndrift::Integer=1,
            ndown::Integer=0,
            PAcutoff::Union{AbstractFloat,Nothing}=nothing)
    drift = par[1:ndrift]
    down = vcat(0.0,par[ndrift+1:ndrift+ndown])
    mfrac = isnothing(mf) ? par[ndrift+ndown+1] : log(mf)
    adrift = isnothing(PAcutoff) ? drift : par[end-ndrift+1:end]
    out = 0.0
    for (refmat,dat) in dats
        (x0,y0,y1) = anchors[refmat]
        Pm,Dm,dm,O,ft,FT,mf,bPt,bDt,bdt =
            LLprep(bP,bD,bd,dat,channels,mfrac,drift,down;
                   PAcutoff=PAcutoff,adrift=adrift)
        out += LL(Pm,Dm,dm,O,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    end
    return out
end
function LL(Pm::AbstractVector,
            Dm::AbstractVector,
            dm::AbstractVector,
            O::Matrix,
            x0::AbstractFloat,
            y0::AbstractFloat,
            y1::AbstractFloat,
            ft::AbstractVector,
            FT::AbstractVector,
            mf::AbstractFloat,
            bPt::AbstractVector,
            bDt::AbstractVector,
            bdt::AbstractVector)
    pred = predict(Pm,Dm,dm,O,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    P = pred[:,"P"]
    D = pred[:,"D"]
    d = pred[:,"d"]
    return sum(@. (((-(FT*P*ft*(y1-y0))/x0)-D*mf*y0+dm)*(O[3,3]*((-(FT*P*ft*(y1-y0))/x0)-D*mf*y0+dm)+O[3,1]*(Pm-P)+(Dm-D)*O[3,2])+(Dm-D)*(O[2,3]*((-(FT*P*ft*(y1-y0))/x0)-D*mf*y0+dm)+O[2,1]*(Pm-P)+(Dm-D)*O[2,2])+(Pm-P)*(O[1,3]*((-(FT*P*ft*(y1-y0))/x0)-D*mf*y0+dm)+O[1,1]*(Pm-P)+(Dm-D)*O[1,2]))/2 )
end
# glass
function LL(par::AbstractVector,
            bD::AbstractVector,
            bd::AbstractVector,
            dats::AbstractDict,
            channels::AbstractDict,
            anchors::AbstractDict)
    mf = exp(par[1])
    out = 0.0
    for (refmat,dat) in dats
        y0 = anchors[refmat]
        Dm,dm,O,bDt,bdt = LLprep(bD,bd,dat,channels)
        out += LL(Dm,dm,O,y0,mf,bDt,bdt)
    end
    return out
end
function LL(Dm::AbstractVector,
            dm::AbstractVector,
            O::Matrix,
            y0::AbstractFloat,
            mf::AbstractFloat,
            bDt::AbstractVector,
            bdt::AbstractVector)
    pred = predict(Dm,dm,O,y0,mf,bDt,bdt)
    D = pred[:,"D"]
    d = pred[:,"d"]
    return sum(@. ((dm-D*mf*y0)*(O[2,2]*(dm-D*mf*y0)+(Dm-D)*O[2,1])+(Dm-D)*(O[1,2]*(dm-D*mf*y0)+(Dm-D)*O[1,1]))/2 )
end
export LL

# minerals
function LLprep(bP::AbstractVector,
                bD::AbstractVector,
                bd::AbstractVector,
                dat::AbstractDataFrame,
                channels::AbstractDict,
                mfrac::AbstractFloat,
                drift::AbstractVector,
                down::AbstractVector;
                PAcutoff::Union{AbstractFloat,Nothing}=nothing,
                adrift::AbstractVector=drift)
    t = dat.t
    T = dat.T
    Pm = dat[:,channels["P"]]
    Dm = dat[:,channels["D"]]
    dm = dat[:,channels["d"]]
    O = O_timeseries(Pm,Dm,dm)
    ft = get_drift(Pm,t,drift;
                   PAcutoff=PAcutoff,adrift=adrift)
    FT = polyFac(down,T)
    mf = exp(mfrac)
    bPt = polyVal(bP,t)
    bDt = polyVal(bD,t)
    bdt = polyVal(bd,t)
    return Pm,Dm,dm,O,ft,FT,mf,bPt,bDt,bdt
end
# glass
function LLprep(bD::AbstractVector,
                bd::AbstractVector,
                dat::AbstractDataFrame,
                channels::AbstractDict)
    t = dat.t
    Dm = dat[:,channels["D"]]
    dm = dat[:,channels["d"]]
    O = O_timeseries(Dm,dm)
    bDt = polyVal(bD,t)
    bdt = polyVal(bd,t)
    return Dm,dm,O,bDt,bdt
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
    Pm,Dm,dm,O,ft,FT,mf,bPt,bDt,bdt =
        LLprep(bP,bD,bd,dat,channels,
               pars.mfrac,pars.drift,pars.down;
               PAcutoff=pars.PAcutoff,adrift=pars.adrift)
    return predict(Pm,Dm,dm,O,
                   anchor.x0,anchor.y0,anchor.y1,
                   ft,FT,mf,bPt,bDt,bdt)
end
function predict(Pm::AbstractVector,
                 Dm::AbstractVector,
                 dm::AbstractVector,
                 O::Matrix,
                 x0::AbstractFloat,
                 y0::AbstractFloat,
                 y1::AbstractFloat,
                 ft::AbstractVector,
                 FT::AbstractVector,
                 mf::AbstractFloat,
                 bPt::AbstractVector,
                 bDt::AbstractVector,
                 bdt::AbstractVector)
    P = getP(Pm,Dm,dm,O,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
    D = getD(Pm,Dm,dm,O,x0,y0,y1,ft,FT,mf,bPt,bDt,bdt)
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
    df = @. D*y0*mf + P*ft*FT*(y1-y0)/x0 + bdt
    return DataFrame(P=Pf,D=Df,d=df)
end
# glass
function predict(dat::AbstractDataFrame,
                 pars::NamedTuple,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 y0::AbstractFloat;
                 debug::Bool=false)
    mf = exp(pars.mfrac)
    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    Dm,dm,O,bDt,bdt = LLprep(bD,bd,dat,channels)
    return predict(Dm,dm,O,y0,mf,bDt,bdt)
end
function predict(Dm::AbstractVector,
                 dm::AbstractVector,
                 O::Matrix,
                 y0::AbstractFloat,
                 mf::AbstractFloat,
                 bDt::AbstractVector,
                 bdt::AbstractVector)
    D = getD(Dm,dm,O,y0,mf,bDt,bdt)
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
# blank
function predict(samp::Sample,
                 blank::AbstractDataFrame;
                 debug::Bool=false)
    dat = windowData(samp;blank=true)
    return polyVal(blank,dat.t)
end
export predict

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
