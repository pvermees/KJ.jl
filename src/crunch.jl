"""
getP(x0::Real,
     y0::Real,
     y1::Real;
     pmb::AbstractVector,
     Domb::AbstractVector,
     bomb::AbstractVector,
     pomb::Union{AbstractVector,Real}=0.0,
     somb::Union{AbstractVector,Real}=0.0,
     vp::Real,
     vD::Real,
     vb::Real,
     spDo::Real,
     spbo::Real,
     sDobo::Real,
     ft::AbstractVector,
     FT::AbstractVector,
     bd::Real=1.0,
     Pp::Real=1.0,
     Ss::Real=1.0)

Estimate the true parent intensity from the measurements for isochron-based standards
"""
function getP(x0::Real,
              y0::Real,
              y1::Real;
              pm::AbstractVector,
              Dom::AbstractVector,
              bom::AbstractVector,
              pomb::Union{AbstractVector,Real}=0.0,
              somb::Union{AbstractVector,Real}=0.0,
              vp::Real,
              vD::Real,
              vb::Real,
              spDo::Real,
              spbo::Real,
              sDobo::Real,
              ft::AbstractVector,
              FT::AbstractVector,
              bd::Real=1.0,
              Pp::Real=1.0,
              Ss::Real=1.0)
    return @. ((((Ss*bd^2*somb+Pp*bd^2*pomb-Domb*bd^2)*vp+bd^2*pmb*spDo)*x0*y0+(bd*bomb*vp-bd*pmb*spbo)*x0)*y1+((FT*bd^2*ft*pmb*vDo+(FT*Ss*bd^2*ft*somb+FT*Pp*bd^2*ft*pomb-Domb*FT*bd^2*ft)*spDo)*x0^2+((-(Ss*bd^2*somb)-Pp*bd^2*pomb+Domb*bd^2)*vp-bd^2*pmb*spDo)*x0)*y0^2+(((-(FT*Ss*bd*ft*somb)-FT*Pp*bd*ft*pomb+Domb*FT*bd*ft)*spbo+FT*bd*bomb*ft*spDo-2*FT*bd*ft*pmb*sDobo)*x0^2+(bd*pmb*spbo-bd*bomb*vp)*x0)*y0+(FT*ft*pmb*vbo-FT*bomb*ft*spbo)*x0^2)/(bd^2*vp*y1^2+((2*FT*bd^2*ft*spDo*x0-2*bd^2*vp)*y0-2*FT*bd*ft*spbo*x0)*y1+(FT^2*bd^2*ft^2*vDo*x0^2-2*FT*bd^2*ft*spDo*x0+bd^2*vp)*y0^2+(2*FT*bd*ft*spbo*x0-2*FT^2*bd*ft^2*sDobo*x0^2)*y0+FT^2*ft^2*vbo*x0^2)
end
export getP

"""
getD(x0::Real,
     y0::Real,
     y1::Real;
     pmb::AbstractVector,
     Domb::AbstractVector,
     bomb::AbstractVector,
     pomb::Union{AbstractVector,Number}=0.0,
     somb::Union{AbstractVector,Number}=0.0,
     vp::Real,
     vD::Real,
     vb::Real,
     spDo::Real,
     spbo::Real,
     sDobo::Real,
     ft::AbstractVector,
     FT::AbstractVector,
     bd::Number=1.0,
     Pp::Number=1.0,
     Ss::Number=1.0)

Estimate the true daughter intensity from the measurements for isochron-based standards
"""
function getD(x0::Real,
              y0::Real,
              y1::Real;
              pmb::AbstractVector,
              Domb::AbstractVector,
              bomb::AbstractVector,
              pomb::Union{AbstractVector,Number}=0.0,
              somb::Union{AbstractVector,Number}=0.0,
              vp::Real,
              vD::Real,
              vb::Real,
              spDo::Real,
              spbo::Real,
              sDobo::Real,
              ft::AbstractVector,
              FT::AbstractVector,
              bd::Number=1.0,
              Pp::Number=1.0,
              Ss::Number=1.0)
    return @. -((((Ss*bd^2*somb+Pp*bd^2*pomb-Domb*bd^2)*vp+bd^2*pmb*spDo)*y1^2+(((FT*bd^2*ft*pmb*vDo+(FT*Ss*bd^2*ft*somb+FT*Pp*bd^2*ft*pomb-Domb*FT*bd^2*ft)*spDo)*x0+(-(2*Ss*bd^2*somb)-2*Pp*bd^2*pomb+2*Domb*bd^2)*vp-2*bd^2*pmb*spDo)*y0+((-(2*FT*Ss*bd*ft*somb)-2*FT*Pp*bd*ft*pomb+2*Domb*FT*bd*ft)*spbo-FT*bd*bomb*ft*spDo-FT*bd*ft*pmb*sDobo)*x0)*y1+(((-(FT*Ss*bd^2*ft*somb)-FT*Pp*bd^2*ft*pomb+Domb*FT*bd^2*ft)*spDo-FT*bd^2*ft*pmb*vDo)*x0+(Ss*bd^2*somb+Pp*bd^2*pomb-Domb*bd^2)*vp+bd^2*pmb*spDo)*y0^2+((-(FT^2*bd*bomb*ft^2*vDo)-FT^2*Ss*bd*ft^2*sDobo*somb+(Domb*FT^2*bd*ft^2-FT^2*Pp*bd*ft^2*pomb)*sDobo)*x0^2+((2*FT*Ss*bd*ft*somb+2*FT*Pp*bd*ft*pomb-2*Domb*FT*bd*ft)*spbo+FT*bd*bomb*ft*spDo+FT*bd*ft*pmb*sDobo)*x0)*y0+((FT^2*Ss*ft^2*somb+FT^2*Pp*ft^2*pomb-Domb*FT^2*ft^2)*vbo+FT^2*bomb*ft^2*sDobo)*x0^2)/(bd^2*vp*y1^2+((2*FT*bd^2*ft*spDo*x0-2*bd^2*vp)*y0-2*FT*bd*ft*spbo*x0)*y1+(FT^2*bd^2*ft^2*vDo*x0^2-2*FT*bd^2*ft*spDo*x0+bd^2*vp)*y0^2+(2*FT*bd*ft*spbo*x0-2*FT^2*bd*ft^2*sDobo*x0^2)*y0+FT^2*ft^2*vbo*x0^2))
end

"""
getD(x::Real,
     y::Real;
     pmb::AbstractVector,
     Domb::AbstractVector,
     bomb::AbstractVector,
     pomb::Union{AbstractVector,Number}=0.0,
     somb::Union{AbstractVector,Number}=0.0,
     vp::Real,
     vD::Real,
     vb::Real,
     spDo::Real,
     spbo::Real,
     sDobo::Real,
     ft::AbstractVector,
     FT::AbstractVector,
     bd::Number=1.0,
     Pp::Number=1.0,
     Ss::Number=1.0)

Estimate the true daughter intensity from the measurements for homogeneous standards
"""
function getD(x::Real,
              y::Real;
              pmb::AbstractVector,
              Domb::AbstractVector,
              bomb::AbstractVector,
              pomb::Union{AbstractVector,Number}=0.0,
              somb::Union{AbstractVector,Number}=0.0,
              vp::Real,
              vD::Real,
              vb::Real,
              spDo::Real,
              spbo::Real,
              sDobo::Real,
              ft::AbstractVector,
              FT::AbstractVector,
              bd::Number=1.0,
              Pp::Number=1.0,
              Ss::Number=1.0)
    return @. (((bd*bomb*vDo+Ss*bd*sDobo*somb+(Pp*bd*pomb-Domb*bd)*sDobo)*vp-bd*pmb*spbo*vDo+(-(Ss*bd*somb)-Pp*bd*pomb+Domb*bd)*spDo*spbo-bd*bomb*spDo^2+bd*pmb*sDobo*spDo)*y+((FT*ft*pmb*vDo+(FT*Ss*ft*somb+FT*Pp*ft*pomb-Domb*FT*ft)*spDo)*vbo-FT*bomb*ft*spbo*vDo+((Domb*FT*ft-FT*Pp*ft*pomb)*sDobo-FT*Ss*ft*sDobo*somb)*spbo+FT*bomb*ft*sDobo*spDo-FT*ft*pmb*sDobo^2)*x+((-(Ss*somb)-Pp*pomb+Domb)*vbo-bomb*sDobo)*vp-pmb*spDo*vbo+(Ss*somb+Pp*pomb-Domb)*spbo^2+(bomb*spDo+pmb*sDobo)*spbo)/((bd^2*vDo*vp-bd^2*spDo^2)*y^2+((2*FT*bd*ft*sDobo*spDo-2*FT*bd*ft*spbo*vDo)*x-2*bd*sDobo*vp+2*bd*spDo*spbo)*y+(FT^2*ft^2*vDo*vbo-FT^2*ft^2*sDobo^2)*x^2+(2*FT*ft*sDobo*spbo-2*FT*ft*spDo*vbo)*x+vbo*vp-spbo^2)
end

"""
SS(x0::Real,
   y0::Real,
   y1::Real;
   pmb::AbstractVector,
   Domb::AbstractVector,
   bomb::AbstractVector,
   pomb::Union{AbstractVector,Number}=0.0,
   somb::Union{AbstractVector,Number}=0.0,
   vp::Real,
   vD::Real,
   vb::Real,
   spDo::Real,
   spbo::Real,
   sDobo::Real,
   ft::AbstractVector,
   FT::AbstractVector,
   bd::Number=1.0,
   Pp::Number=1.0,
   Ss::Number=1.0)

Sum of squares for isochron-based standards
"""
function SS(x0::Real,
            y0::Real,
            y1::Real;
            pmb::AbstractVector,
            Domb::AbstractVector,
            bomb::AbstractVector,
            pomb::Union{AbstractVector,Number}=0.0,
            somb::Union{AbstractVector,Number}=0.0,
            vp::Real,
            vD::Real,
            vb::Real,
            spDo::Real,
            spbo::Real,
            sDobo::Real,
            ft::AbstractVector,
            FT::AbstractVector,
            bd::Number=1.0,
            Pp::Number=1.0,
            Ss::Number=1.0)
    P = getP(x0,y0,y1;
             pmb=pmb,Domb=Domb,bomb=bomb,pomb=pomb,somb=somb,
             vp=vp,vD=vD,vb=vb,spDo=spDo,spbo=spbo,sDobo=sDobo,
             ft=ft,FT=FT,bd=bd,Pp=Pp,Ss=Ss)
    D = getD(x0,y0,y1;
             pmb=pmb,Domb=Domb,bomb=bomb,pomb=pomb,somb=somb,
             vp=vp,vD=vD,vb=vb,spDo=spDo,spbo=spbo,sDobo=sDobo,
             ft=ft,FT=FT,bd=bd,Pp=Pp,Ss=Ss)
    maha = @. (bomb-bd*((Po*(y1-y0))/x0+Do*y0))*(((vDo*vp-spDo^2)*(bomb-bd*((Po*(y1-y0))/x0+Do*y0)))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo))+((-(Ss*somb)-Pp*pomb+Domb-Do)*(spDo*spbo-sDobo*vp))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo))+((pmb-FT*Po*ft)*(sDobo*spDo-spbo*vDo))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo)))+(-(Ss*somb)-Pp*pomb+Domb-Do)*(((spDo*spbo-sDobo*vp)*(bomb-bd*((Po*(y1-y0))/x0+Do*y0)))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo))+((-(Ss*somb)-Pp*pomb+Domb-Do)*(vbo*vp-spbo^2))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo))+((pmb-FT*Po*ft)*(sDobo*spbo-spDo*vbo))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo)))+(pmb-FT*Po*ft)*(((sDobo*spDo-spbo*vDo)*(bomb-bd*((Po*(y1-y0))/x0+Do*y0)))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo))+((pmb-FT*Po*ft)*(vDo*vbo-sDobo^2))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo))+((-(Ss*somb)-Pp*pomb+Domb-Do)*(sDobo*spbo-spDo*vbo))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo)))
    return sum(@. maha )
end

"""
SS(x::Real,
   y::Real;
   pmb::AbstractVector,
   Domb::AbstractVector,
   bomb::AbstractVector,
   pomb::Union{AbstractVector,Number}=0.0,
   somb::Union{AbstractVector,Number}=0.0,
   vp::Real,
   vD::Real,
   vb::Real,
   spDo::Real,
   spbo::Real,
   sDobo::Real,
   ft::AbstractVector,
   FT::AbstractVector,
   bd::Number=1.0,
   Pp::Number=1.0,
   Ss::Number=1.0)

Sum of squares for homogeneous standards
"""
function SS(x::Real,
            y::Real;
            pmb::AbstractVector,
            Domb::AbstractVector,
            bomb::AbstractVector,
            pomb::Union{AbstractVector,Number}=0.0,
            somb::Union{AbstractVector,Number}=0.0,
            vp::Real,
            vD::Real,
            vb::Real,
            spDo::Real,
            spbo::Real,
            sDobo::Real,
            ft::AbstractVector,
            FT::AbstractVector,
            bd::Number=1.0,
            Pp::Number=1.0,
            Ss::Number=1.0)
    D = getD(x,y;
             pmb=pmb,Domb=Domb,bomb=bomb,pomb=pomb,somb=somb,
             vp=vp,vD=vD,vb=vb,spDo=spDo,spbo=spbo,sDobo=sDobo,
             ft=ft,FT=FT,bd=bd,Pp=Pp,Ss=Ss)
    maha = @. (bomb-Do*bd*y)*(((vDo*vp-spDo^2)*(bomb-Do*bd*y))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo))+((sDobo*spDo-spbo*vDo)*(pmb-Do*FT*ft*x))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo))+((-(Ss*somb)-Pp*pomb+Domb-Do)*(spDo*spbo-sDobo*vp))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo)))+(-(Ss*somb)-Pp*pomb+Domb-Do)*(((spDo*spbo-sDobo*vp)*(bomb-Do*bd*y))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo))+((sDobo*spbo-spDo*vbo)*(pmb-Do*FT*ft*x))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo))+((-(Ss*somb)-Pp*pomb+Domb-Do)*(vbo*vp-spbo^2))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo)))+(pmb-Do*FT*ft*x)*(((sDobo*spDo-spbo*vDo)*(bomb-Do*bd*y))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo))+((vDo*vbo-sDobo^2)*(pmb-Do*FT*ft*x))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo))+((-(Ss*somb)-Pp*pomb+Domb-Do)*(sDobo*spbo-spDo*vbo))/((vDo*vbo-sDobo^2)*vp+spDo*(sDobo*spbo-spDo*vbo)+spbo*(sDobo*spDo-spbo*vDo)))
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
        dat = swinData(samp)
        sig = getSignals(dat)
        covmat = df2cov(sig)
        return predict(dat,covmat,method,fit;kw...)
    else
        KJerror("notStandard")
    end
end

"""
predict(dat::AbstractDataFrame,
        covmat::Matrix,
        method::KJmethod,
        fit::KJfit)
"""
function predict(dat::AbstractDataFrame,
                 covmat::Matrix,
                 method::KJmethod,
                 fit::KJfit;
                 kw...)
    t = dat.t
    T = dat.T
    ch = getChannels(method)
    sig = getSignals(dat)
    iP = columnindex(sig,ch.P)
    iD = columnindex(sig,ch.D)
    id = columnindex(sig,ch.d)
    Pm = sig[:,iP]
    Dm = sig[:,iD]
    dm = sig[:,id]
    vP = covmat[iP,iP]
    vD = covmat[iD,iD]
    vd = covmat[id,id]
    sPD = covmat[iP,iD]
    sPd = covmat[iP,id]
    sDd = covmat[iD,id]
    ft = get_drift(Pm,t,fit)
    FT = polyFac(fit.down,T)
    ions = getIons(method)
    proxies = getProxies(method)
    bd = iratio(ions.d,proxies.d)
    Ss = iratio(ions.S,proxies.S)
    blank = fit.blank
    bPt = polyVal(blank[:,ch.P],t)
    bDt = polyVal(blank[:,ch.D],t)
    bdt = polyVal(blank[:,ch.d],t)
    if is_isochron_anchor(anchor)
        return predict(anchor.x0,
                       anchor.y0,
                       anchor.y1;
                       pm=Pm,Dom=Dm,bom=dm,
                       vp=vP,vD=vD,vd=vd,
                       spDo=sPD,spbo=sPd,sDobo=sDd,
                       ft=ft,FT=FT,bd=bd,Ss=Ss,
                       bpt=bPt,bDt=bDt,bbt=bdt)
    elseif is_point_anchor(anchor)
        return predict(anchor.x,
                       anchor.y;
                       pm=Pm,Dom=Dm,bom=dm,
                       vp=vP,vD=vD,vd=vd,
                       spDo=sPD,spbo=sPd,sDobo=sDd,
                       ft=ft,FT=FT,bd=bd,Ss=Ss,
                       bpt=bPt,bDt=bDt,bbt=bdt)
    else
        error("Invalid anchor")
    end
end

"""
predict(x0::Real,
        y0::Real,
        y1::Real=0.0;
        pm::AbstractVector,
        Dom::AbstractVector,
        bom::AbstractVector,
        pomb::Union{AbstractVector,Real}=0.0,
        somb::Union{AbstractVector,Real}=0.0,
        vp::Real,
        vD::Real,
        vb::Real,
        spDo::Real,
        spbo::Real,
        sDobo::Real,
        ft::AbstractVector,
        FT::AbstractVector,
        bd::Real=1.0,
        Pp::Real=1.0,
        Ss::Real=1.0,
        bpt::AbstractVector,
        bDot::AbstractVector,
        bbot::AbstractVector)

For isochron-based standards
"""
function predict(x0::Real,
                 y0::Real,
                 y1::Real=0.0;
                 pm::AbstractVector,
                 Dom::AbstractVector,
                 bom::AbstractVector,
                 pomb::Union{AbstractVector,Real}=0.0,
                 somb::Union{AbstractVector,Real}=0.0,
                 vp::Real,
                 vD::Real,
                 vb::Real,
                 spDo::Real,
                 spbo::Real,
                 sDobo::Real,
                 ft::AbstractVector,
                 FT::AbstractVector,
                 bd::Real=1.0,
                 Pp::Real=1.0,
                 Ss::Real=1.0,
                 bpt::AbstractVector,
                 bDot::AbstractVector,
                 bbot::AbstractVector)
    Po = getP(x0,y0,y1;
              pm=pm,Dom=Dom,bom=bom,pomb=pomb,somb=somb,
              vp=vp,vD=vD,vb=vb,spDo=spDo,spbo=spbo,sDobo=sDobo,
              ft=ft,FT=FT,bd=bd,Pp=Pp,Ss=Ss)
    Do = getD(x0,y0,y1;
              pm=pm,Dom=Dom,bom=bom,pomb=pomb,somb=somb,
              vp=vp,vD=vD,vb=vb,spDo=spDo,spbo=spbo,sDobo=sDobo,
              ft=ft,FT=FT,bd=bd,Pp=Pp,Ss=Ss)
    Pf = @. Po*ft*FT + bpt
    Dof = @. Do + pomb*Pp + somb*Ss + bDot
    bof = @. (Do*y0 + Po*(y1-y0)/x0)*bd + bbot
    return DataFrame(P=Pf,D=Dof,d=bof)
end

"""
predict(x::Real,
        y::Real;
        pm::AbstractVector,
        Dom::AbstractVector,
        bom::AbstractVector,
        pomb::Union{AbstractVector,Real}=0.0,
        somb::Union{AbstractVector,Real}=0.0,
        vp::Real,
        vD::Real,
        vb::Real,
        spDo::Real,
        spbo::Real,
        sDobo::Real,
        ft::AbstractVector,
        FT::AbstractVector,
        bd::Real=1.0,
        Pp::Real=1.0,
        Ss::Real=1.0,
        bPt::AbstractVector,
        bDot::AbstractVector,
        bbot::AbstractVector)

For homogeneous standards
"""
function predict(x::Real,
                 y::Real;
                 pm::AbstractVector,
                 Dom::AbstractVector,
                 bom::AbstractVector,
                 pomb::Union{AbstractVector,Real}=0.0,
                 somb::Union{AbstractVector,Real}=0.0,
                 vp::Real,
                 vD::Real,
                 vb::Real,
                 spDo::Real,
                 spbo::Real,
                 sDobo::Real,
                 ft::AbstractVector,
                 FT::AbstractVector,
                 bd::Real=1.0,
                 Pp::Real=1.0,
                 Ss::Real=1.0,
                 bPt::AbstractVector,
                 bDot::AbstractVector,
                 bbot::AbstractVector)
    Do = getD(x,y;
              pm=pm,Dom=Dom,bom=bom,pomb=pomb,somb=somb,
              vp=vp,vD=vD,vb=vb,spDo=spDo,spbo=spbo,sDobo=sDobo,
              ft=ft,FT=FT,bd=bd,Pp=Pp,Ss=Ss)
    Pf = @. Do*x*ft*FT + bPt
    Dof = @. Do + pomb*Pp + somb*Ss + bDot
    bof = @. Do*y*bd + bdot
    return DataFrame(P=Pf,D=Dof,d=bof)
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
