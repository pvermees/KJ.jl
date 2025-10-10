# two-step isotope fractionation
function fractionation(run::Vector{Sample},
                       method::AbstractString,
                       blank::AbstractDataFrame,
                       channels::AbstractDict,
                       standards::AbstractDict,
                       glass::AbstractDict;
                       ndrift::Integer=1,
                       ndown::Integer=0,
                       PAcutoff=nothing,
                       verbose::Bool=false)
    mf = fractionation(run,method,blank,channels,glass;verbose=verbose)
    return fractionation(run,method,blank,channels,standards,mf;
                         ndrift=ndrift,ndown=ndown,
                         PAcutoff=PAcutoff,verbose=verbose)
end
# one-step isotope fractionation using mineral standards
function fractionation(run::Vector{Sample},
                       method::AbstractString,
                       blank::AbstractDataFrame,
                       channels::AbstractDict,
                       standards::AbstractDict,
                       mf::Union{Real,Nothing};
                       ndrift::Integer=1,
                       ndown::Integer=0,
                       PAcutoff=nothing,
                       verbose::Bool=false)
    
    anchors = getStandardAnchors(method,standards)
    
    if ndrift<1 KJerror("ndriftzero") end

    dats = Dict()
    covs = Dict()
    for (refmat,anchor) in anchors
        dats[refmat], covs[refmat] = pool(run;
                                          signal=true,
                                          group=refmat,
                                          include_covmats=true)
    end

    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    bP = blank[:,channels["P"]]

    init = fill(0.0,ndrift)
    if (ndown>0) init = vcat(init,fill(0.0,ndown)) end
    if isnothing(mf) init = vcat(init,0.0) end
    if !isnothing(PAcutoff) init = vcat(init,fill(0.0,ndrift)) end
    
    return SSfit(init,bP,bD,bd,dats,covs,channels,anchors,mf;
                 ndrift=ndrift,ndown=ndown,
                 PAcutoff=PAcutoff,verbose=verbose)

end
# isotopic mass fractionation using glass
function fractionation(run::Vector{Sample},
                       method::AbstractString,
                       blank::AbstractDataFrame,
                       channels::AbstractDict,
                       glass::AbstractDict;
                       verbose::Bool=false)
    
    anchors = getGlassAnchors(method,glass)

    dats = Dict()
    covs = Dict()
    for (refmat,anchor) in anchors
        dats[refmat], covs[refmat] = pool(run;signal=true,group=refmat,
                                          include_covmats=true)
    end

    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]

    return SSfit([0.0],bD,bd,dats,covs,channels,anchors;verbose=verbose)
    
end
# for concentration measurements:
function fractionation(run::Vector{Sample},
                       blank::AbstractDataFrame,
                       internal::Tuple,
                       glass::AbstractDict)
    elements = channels2elements(run)
    return fractionation(run,blank,elements,internal,glass)
end
function fractionation(run::Vector{Sample},
                       blank::AbstractDataFrame,
                       elements::AbstractDataFrame,
                       internal::Tuple,
                       glass::AbstractDict)
    ne = size(elements,2)
    num = den = fill(0.0,ne-1)
    for (SRM,name)  in glass
        dat = pool(run;signal=true,group=SRM)
        concs = elements2concs(elements,SRM)
        bt = polyVal(blank,dat.t)
        sig = getSignals(dat)
        (nr,nc) = size(sig)
        Xm = sig[:,Not(internal[1])]
        Sm = sig[:,internal[1]]
        bXt = bt[:,Not(internal[1])]
        bSt = bt[:,internal[1]]
        S = Sm.-bSt
        R = collect((concs[:,Not(internal[1])]./concs[:,internal[1]])[1,:])
        num += sum.(eachcol(R'.*(Xm.-bXt).*S))
        den += sum.(eachcol((R'.*S).^2))
    end
    return num./den
end
export fractionation

# minerals
function SSfit(init::AbstractVector,
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
               PAcutoff=nothing,
               verbose::Bool=false)

    objective = (par) -> SS(par,bP,bD,bd,dats,covs,channels,anchors,mf;
                            ndrift=ndrift,ndown=ndown,PAcutoff=PAcutoff)

    fit = Optim.optimize(objective,init)
    pars = Optim.minimizer(fit)
    if verbose
        println("Drift and downhole fractionation correction:\n")
        println(fit)
    else
        if fit.stopped_by.time_limit
            @warn "Reached the maximum number of iterations before achieving " *
                "convergence. Reduce the order of the polynomials or fix the " *
                "mass fractionation and try again."
        end
        if hasproperty(fit.stopped_by,:ls_failed) && fit.stopped_by.ls_failed
            @warn "Least squares algorithm did not converge."
        end
    end
    drift = pars[1:ndrift]
    down = vcat(0.0,pars[ndrift+1:ndrift+ndown])
    mfrac = isnothing(mf) ? pars[ndrift+ndown+1] : log(mf)
    adrift = isnothing(PAcutoff) ? drift : pars[end-ndrift+1:end]
    
    return (drift=drift,down=down,mfrac=mfrac,
            PAcutoff=PAcutoff,adrift=adrift)

end
# glass
function SSfit(init::AbstractVector,
               bD::AbstractVector,
               bd::AbstractVector,
               dats::AbstractDict,
               vars::AbstractDict,
               channels::AbstractDict,
               anchors::AbstractDict;
               verbose::Bool=false)

    objective = (par) -> SS(par,bD,bd,dats,vars,channels,anchors)
    
    fit = Optim.optimize(objective,init)
    pars = Optim.minimizer(fit)

    if verbose
        println("Mass fractionation correction:\n")
        println(fit)
    end
    
    mfrac = Optim.minimizer(fit)[1]
        
    return exp(mfrac)
    
end
