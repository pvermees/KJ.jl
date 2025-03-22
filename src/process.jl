# ratios:
function process!(run::Vector{Sample},
                  method::AbstractString,
                  channels::AbstractDict,
                  standards::AbstractDict,
                  glass::AbstractDict;
                  nblank::Integer=2,ndrift::Integer=1,ndown::Integer=1,
                  PAcutoff=nothing,verbose::Bool=false)
    blank = fitBlanks(run;nblank=nblank)
    setGroup!(run,glass)
    setGroup!(run,standards)
    fit = fractionation(run,method,blank,channels,standards,glass;
                        ndrift=ndrift,ndown=ndown,
                        PAcutoff=PAcutoff,verbose=verbose)
    return blank, fit
end
# concentrations:
function process!(run::Vector{Sample},
                  internal::Tuple,
                  glass::AbstractDict;
                  nblank::Integer=2)
    blank = fitBlanks(run;nblank=nblank)
    setGroup!(run,glass)
    fit = fractionation(run,blank,internal,glass)
    return blank, fit
end
export process!

function fitBlanks(run::Vector{Sample};
                   nblank=2)
    blk = pool(run;blank=true)
    channels = getChannels(run)
    nc = length(channels)
    bpar = DataFrame(zeros(nblank,nc),channels)
    for channel in channels
        bpar[:,channel] = polyFit(blk.t,blk[:,channel],nblank)
    end
    return bpar
end
export fitBlanks

function atomic(samp::Sample,
                channels::AbstractDict,
                blank::AbstractDataFrame,
                pars::NamedTuple)
    dat = windowData(samp;signal=true)
    var = dat2var(dat,collect(values(channels)))
    Pm,Dm,dm,vP,vD,vd,ft,FT,mf,bPt,bDt,bdt =
        SSprep(blank[:,channels["P"]],
               blank[:,channels["D"]],
               blank[:,channels["d"]],
               dat,var,
               channels,
               pars.mfrac,pars.drift,pars.down;
               PAcutoff=pars.PAcutoff,
               adrift=pars.adrift)
    Phat = @. (Pm-bPt)*ft*FT
    Dhat = @. (Dm-bDt)*mf
    dhat = @. (dm-bdt)
    return Phat, Dhat, dhat
end
export atomic

function concentrations(samp::Sample,
                        blank::AbstractDataFrame,
                        pars::AbstractVector,
                        internal::Tuple)
    elements = channels2elements(samp)
    return concentrations(samp,elements,blank,pars,internal)
end
function concentrations(samp::Sample,
                        elements::AbstractDataFrame,
                        blank::AbstractDataFrame,
                        pars::AbstractVector,
                        internal::Tuple)
    dat = windowData(samp;signal=true)
    sig = getSignals(dat)
    out = copy(sig)
    bt = polyVal(blank,dat.t)
    bXt = bt[:,Not(internal[1])]
    bSt = bt[:,internal[1]]
    Xm = sig[:,Not(internal[1])]
    Sm = sig[:,internal[1]]
    out[!,internal[1]] .= internal[2]
    num = @. (Xm-bXt)*internal[2]
    den = @. pars'*(Sm-bSt)
    out[!,Not(internal[1])] .= num./den
    elementnames = collect(elements[1,:])
    channelnames = names(sig)
    nms = "ppm[" .* elementnames .* "] from " .* channelnames
    rename!(out,Symbol.(nms))
    return out
end
function concentrations(run::Vector{Sample},
                        blank::AbstractDataFrame,
                        pars::AbstractVector,
                        internal::Tuple)
    elements = channels2elements(run)
    return concentrations(run,elements,blank,pars,internal)
end
function concentrations(run::Vector{Sample},
                        elements::AbstractDataFrame,
                        blank::AbstractDataFrame,
                        pars::AbstractVector,
                        internal::Tuple)
    nr = length(run)
    nc = 2*size(elements,2)
    mat = fill(0.0,nr,nc)
    conc = nothing
    for i in eachindex(run)
        samp = run[i]
        conc = concentrations(samp,elements,blank,pars,internal)
        mu = Statistics.mean.(eachcol(conc))
        sigma = Statistics.std.(eachcol(conc))
        mat[i,1:2:nc-1] .= mu
        mat[i,2:2:nc] .= sigma
    end
    nms = fill("",nc)
    nms[1:2:nc-1] .= names(conc)
    nms[2:2:nc] .= "s[" .* names(conc) .* "]"
    out = hcat(DataFrame(sample=getSnames(run)),DataFrame(mat,Symbol.(nms)))
    return out
end
export concentrations
