function getP(a::IsochronAnchor,
              ft::AbstractVector,
              hT::AbstractVector;
              pmb::AbstractVector,
              Dmb::AbstractVector,
              bmb::AbstractVector,
              vp::AbstractFloat,
              vD::AbstractFloat,
              vb::AbstractFloat,
              spD::AbstractFloat,
              spb::AbstractFloat,
              sDb::AbstractFloat,
              mf::AbstractVector,
              bd::AbstractFloat,
              Ip::AbstractVector,
              ID::AbstractVector,
              Ib::AbstractVector,
              other...)
    x0,y0,y1 = unpack(a)
    return @. ((((ID-Dmb)*bd^2*mf^2*vp+(bd^2*mf^2*pmb-Ip*bd^2*mf^2)*spD)*x0*y0+((bd*bmb-Ib*bd)*mf*vp+(Ip*bd*mf-bd*mf*pmb)*spb)*x0)*y1+(((bd^2*ft*hT*mf^2*pmb-Ip*bd^2*ft*hT*mf^2)*vD+(ID-Dmb)*bd^2*ft*hT*mf^2*spD)*x0^2+((Dmb-ID)*bd^2*mf^2*vp+(Ip*bd^2*mf^2-bd^2*mf^2*pmb)*spD)*x0)*y0^2+(((Dmb-ID)*bd*ft*hT*mf*spb+(bd*bmb-Ib*bd)*ft*hT*mf*spD+(2*Ip*bd*ft*hT*mf-2*bd*ft*hT*mf*pmb)*sDb)*x0^2+((Ib*bd-bd*bmb)*mf*vp+(bd*mf*pmb-Ip*bd*mf)*spb)*x0)*y0+((ft*hT*pmb-Ip*ft*hT)*vb+(Ib-bmb)*ft*hT*spb)*x0^2)/(bd^2*mf^2*vp*y1^2+((2*bd^2*ft*hT*mf^2*spD*x0-2*bd^2*mf^2*vp)*y0-2*bd*ft*hT*mf*spb*x0)*y1+(bd^2*ft^2*hT^2*mf^2*vD*x0^2-2*bd^2*ft*hT*mf^2*spD*x0+bd^2*mf^2*vp)*y0^2+(2*bd*ft*hT*mf*spb*x0-2*bd*ft^2*hT^2*mf*sDb*x0^2)*y0+ft^2*hT^2*vb*x0^2)
end
export getP

function getD(a::IsochronAnchor,
              ft::AbstractVector,
              hT::AbstractVector;
              pmb::AbstractVector,
              Dmb::AbstractVector,
              bmb::AbstractVector,
              vp::AbstractFloat,
              vD::AbstractFloat,
              vb::AbstractFloat,
              spD::AbstractFloat,
              spb::AbstractFloat,
              sDb::AbstractFloat,
              mf::AbstractVector,
              bd::AbstractFloat,
              Ip::AbstractVector,
              ID::AbstractVector,
              Ib::AbstractVector,
              other...)
    x0,y0,y1 = unpack(a)
    return @. -((((ID-Dmb)*bd^2*mf^2*vp+(bd^2*mf^2*pmb-Ip*bd^2*mf^2)*spD)*y1^2+((((bd^2*ft*hT*mf^2*pmb-Ip*bd^2*ft*hT*mf^2)*vD+(ID-Dmb)*bd^2*ft*hT*mf^2*spD)*x0+(2*Dmb-2*ID)*bd^2*mf^2*vp+(2*Ip*bd^2*mf^2-2*bd^2*mf^2*pmb)*spD)*y0+((2*Dmb-2*ID)*bd*ft*hT*mf*spb+(Ib*bd-bd*bmb)*ft*hT*mf*spD+(Ip*bd*ft*hT*mf-bd*ft*hT*mf*pmb)*sDb)*x0)*y1+(((Ip*bd^2*ft*hT*mf^2-bd^2*ft*hT*mf^2*pmb)*vD+(Dmb-ID)*bd^2*ft*hT*mf^2*spD)*x0+(ID-Dmb)*bd^2*mf^2*vp+(bd^2*mf^2*pmb-Ip*bd^2*mf^2)*spD)*y0^2+(((Ib*bd-bd*bmb)*ft^2*hT^2*mf*vD+(Dmb-ID)*bd*ft^2*hT^2*mf*sDb)*x0^2+((2*ID-2*Dmb)*bd*ft*hT*mf*spb+(bd*bmb-Ib*bd)*ft*hT*mf*spD+(bd*ft*hT*mf*pmb-Ip*bd*ft*hT*mf)*sDb)*x0)*y0+((ID-Dmb)*ft^2*hT^2*vb+(bmb-Ib)*ft^2*hT^2*sDb)*x0^2)/(bd^2*mf^2*vp*y1^2+((2*bd^2*ft*hT*mf^2*spD*x0-2*bd^2*mf^2*vp)*y0-2*bd*ft*hT*mf*spb*x0)*y1+(bd^2*ft^2*hT^2*mf^2*vD*x0^2-2*bd^2*ft*hT*mf^2*spD*x0+bd^2*mf^2*vp)*y0^2+(2*bd*ft*hT*mf*spb*x0-2*bd*ft^2*hT^2*mf*sDb*x0^2)*y0+ft^2*hT^2*vb*x0^2))
end

function getD(a::PointAnchor,
              ft::AbstractVector,
              hT::AbstractVector;
              pmb::AbstractVector,
              Dmb::AbstractVector,
              bmb::AbstractVector,
              vp::AbstractFloat,
              vD::AbstractFloat,
              vb::AbstractFloat,
              spD::AbstractFloat,
              spb::AbstractFloat,
              sDb::AbstractFloat,
              mf::AbstractVector,
              bd::AbstractFloat,
              Ip::AbstractVector,
              ID::AbstractVector,
              Ib::AbstractVector,
              other...)
    x,y = unpack(a)
    return @. ((((bd*bmb-Ib*bd)*mf*vD+(ID-Dmb)*bd*mf*sDb)*vp+(Ip*bd*mf-bd*mf*pmb)*spb*vD+(Dmb-ID)*bd*mf*spD*spb+(Ib*bd-bd*bmb)*mf*spD^2+(bd*mf*pmb-Ip*bd*mf)*sDb*spD)*y+(((ft*hT*pmb-Ip*ft*hT)*vD+(ID-Dmb)*ft*hT*spD)*vb+(Ib-bmb)*ft*hT*spb*vD+(Dmb-ID)*ft*hT*sDb*spb+(bmb-Ib)*ft*hT*sDb*spD+(Ip*ft*hT-ft*hT*pmb)*sDb^2)*x+((Dmb-ID)*vb+(Ib-bmb)*sDb)*vp+(Ip-pmb)*spD*vb+(ID-Dmb)*spb^2+((bmb-Ib)*spD+(pmb-Ip)*sDb)*spb)/((bd^2*mf^2*vD*vp-bd^2*mf^2*spD^2)*y^2+((2*bd*ft*hT*mf*sDb*spD-2*bd*ft*hT*mf*spb*vD)*x-2*bd*mf*sDb*vp+2*bd*mf*spD*spb)*y+(ft^2*hT^2*vD*vb-ft^2*hT^2*sDb^2)*x^2+(2*ft*hT*sDb*spb-2*ft*hT*spD*vb)*x+vb*vp-spb^2)
end

function getD(mf::AbstractVector,
              y::AbstractFloat,
              bd::AbstractFloat;
              Dmb::AbstractVector,
              bmb::AbstractVector,
              vD::AbstractFloat,
              vb::AbstractFloat,
              sDb::AbstractFloat,
              other...)
    return @. ((bd*bmb*mf*vD-Dmb*bd*mf*sDb)*y+Dmb*vb-bmb*sDb)/(bd^2*mf^2*vD*y^2-2*bd*mf*sDb*y+vb)
end

function mahalanobis(a::IsochronAnchor,
                     ft::AbstractVector,
                     hT::AbstractVector,
                     P::AbstractVector,
                     D::AbstractVector;
                     pmb::AbstractVector,
                     Dmb::AbstractVector,
                     bmb::AbstractVector,
                     vp::AbstractFloat,
                     vD::AbstractFloat,
                     vb::AbstractFloat,
                     spD::AbstractFloat,
                     spb::AbstractFloat,
                     sDb::AbstractFloat,
                     mf::AbstractVector,
                     bd::AbstractFloat,
                     Ip::AbstractVector,
                     ID::AbstractVector,
                     Ib::AbstractVector,
                     other...)
    x0,y0,y1 = unpack(a)
    return @. (-(bd*mf*((P*(y1-y0))/x0+D*y0))+bmb-Ib)*(((vD*vp-spD^2)*(-(bd*mf*((P*(y1-y0))/x0+D*y0))+bmb-Ib))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((-ID+Dmb-D)*(spD*spb-sDb*vp))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((pmb-P*ft*hT-Ip)*(sDb*spD-spb*vD))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(-ID+Dmb-D)*(((spD*spb-sDb*vp)*(-(bd*mf*((P*(y1-y0))/x0+D*y0))+bmb-Ib))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((-ID+Dmb-D)*(vb*vp-spb^2))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((pmb-P*ft*hT-Ip)*(sDb*spb-spD*vb))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(pmb-P*ft*hT-Ip)*(((sDb*spD-spb*vD)*(-(bd*mf*((P*(y1-y0))/x0+D*y0))+bmb-Ib))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((pmb-P*ft*hT-Ip)*(vD*vb-sDb^2))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((-ID+Dmb-D)*(sDb*spb-spD*vb))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))
end

function mahalanobis(a::PointAnchor,
                     ft::AbstractVector,
                     hT::AbstractVector,
                     D::AbstractVector;
                     pmb::AbstractVector,
                     Dmb::AbstractVector,
                     bmb::AbstractVector,
                     vp::AbstractFloat,
                     vD::AbstractFloat,
                     vb::AbstractFloat,
                     spD::AbstractFloat,
                     spb::AbstractFloat,
                     sDb::AbstractFloat,
                     mf::AbstractVector,
                     bd::AbstractFloat,
                     Ip::AbstractVector,
                     ID::AbstractVector,
                     Ib::AbstractVector,
                     other...)
    x,y= unpack(a)
    return @. (-(D*bd*mf*y)+bmb-Ib)*(((vD*vp-spD^2)*(-(D*bd*mf*y)+bmb-Ib))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((sDb*spD-spb*vD)*(-(D*ft*hT*x)+pmb-Ip))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((-ID+Dmb-D)*(spD*spb-sDb*vp))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(-ID+Dmb-D)*(((spD*spb-sDb*vp)*(-(D*bd*mf*y)+bmb-Ib))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((sDb*spb-spD*vb)*(-(D*ft*hT*x)+pmb-Ip))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((-ID+Dmb-D)*(vb*vp-spb^2))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))+(-(D*ft*hT*x)+pmb-Ip)*(((sDb*spD-spb*vD)*(-(D*bd*mf*y)+bmb-Ib))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((vD*vb-sDb^2)*(-(D*ft*hT*x)+pmb-Ip))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD))+((-ID+Dmb-D)*(sDb*spb-spD*vb))/((vD*vb-sDb^2)*vp+spD*(sDb*spb-spD*vb)+spb*(sDb*spD-spb*vD)))
end

function mahalanobis(mf::AbstractVector,
                     y::AbstractFloat,
                     bd::AbstractFloat,
                     D::AbstractVector;
                     Dmb::AbstractVector,
                     bmb::AbstractVector,
                     vD::AbstractFloat,
                     vb::AbstractFloat,
                     sDb::AbstractFloat,
                     other...)
    return @. (bmb-D*bd*mf*y)*((vD*(bmb-D*bd*mf*y))/(vD*vb-sDb^2)-((Dmb-D)*sDb)/(vD*vb-sDb^2))+(Dmb-D)*(((Dmb-D)*vb)/(vD*vb-sDb^2)-(sDb*(bmb-D*bd*mf*y))/(vD*vb-sDb^2))
end

function SS(a::IsochronAnchor,
            ft::AbstractVector,
            hT::AbstractVector;
            cruncher...)
    P = getP(a,ft,hT;cruncher...)
    D = getD(a,ft,hT;cruncher...)
    maha = mahalanobis(a,ft,hT,P,D;cruncher...)
    return sum(@. maha )
end

function SS(a::PointAnchor,
            ft::AbstractVector,
            hT::AbstractVector;
            cruncher...)
    D = getD(a,ft,hT;cruncher...)
    maha = mahalanobis(a,ft,hT,D;cruncher...)
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
            ft, hT = ft_hT(fit,method.PAcutoff;cruncher...)
            out += SS(a,ft,hT;cruncher...)
        end
    end
    if verbose
        println(par,": ",out)
    end
    return out
end

function SS(bias::Bias,
            y::AbstractFloat,
            bd::AbstractFloat;
            cruncher...)
    mf = bias_correction(bias,bias.mass_num,bias.mass_den;cruncher...)
    D = getD(mf,y,bd;cruncher...)
    maha = mahalanobis(mf,y,bd,D;cruncher...)
    return sum(@. maha )
end
function SS(bias::REEBias;
            cruncher...)
    mf = bias_correction(bias;cruncher...)
    D = getD(mf,1.0,1.0;cruncher...)
    maha = mahalanobis(mf,1.0,1.0,D;cruncher...)
    return sum(@. maha )
end

function SS(par::AbstractVector,
            mass_num::Int,
            mass_den::Int,
            bd::AbstractFloat,
            cruncher_groups::AbstractDict;
            verbose::Bool=false)
    out = 0.0
    bias = Bias(mass_num,mass_den,par)
    for cg in values(cruncher_groups)
        for cruncher in cg.crunchers
            out += SS(bias,cg.y,bd;cruncher...)
        end
    end
    if verbose
        println(par,": ",out)
    end
    return out
end
function SS(par::AbstractVector,
            cruncher_groups::AbstractDict;
            verbose::Bool=false)
    out = 0.0
    bias = REEBias(par)
    for cruncher_group in values(cruncher_groups)
        for cruncher in cruncher_group
            out += SS(bias;cruncher...)
        end
    end
    if verbose
        println(par,": ",out)
    end
    return out
end
export SS