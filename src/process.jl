"""
process!

Fits blanks and fractionation effects

# Returns

- `blank`: a dataframe with blank parameters for all specified channels
- `fit`: a tuple or (if method is omitted) vector of fit parameters

# Methods

- `process!(run::Vector{Sample},
            method::AbstractString,
            channels::AbstractDict,
            standards::AbstractDict,
            glass::AbstractDict;
            nblank::Integer=2,ndrift::Integer=1,ndown::Integer=1,
            PAcutoff=nothing,verbose::Bool=false)`
- `process!(run::Vector{Sample},
            internal::Tuple,
            glass::AbstractDict;
            nblank::Integer=2)`

# Arguments

- `run`: the output of `load`
- `method`: either "U-Pb", "Lu-Hf", "Rb-Sr" or "concentrations"
- `channels`: dictionary of the type Dict("P" => "parent", "D" => "daughter", "d" => "sister")
              or a vector of channel names (e.g., the keys of a channels Dict)
- `standards`: dictionary of the type Dict("prefix" => "mineral standard")
- `glass`: dictionary of the type Dict("prefix" => "reference glass")
- `nblank`, `ndrift`, `ndown`: The number of parameters used to fit the blanks,
                               drift and down hole fractionation, respectively
- `PAcutoff`: pulse-analog cutoff
- `verbose`: if `true`, prints the optimisation results to the REPL
- `internal`: a tuple with the name of a channel that is to be used as an internal
              concentration standard, and its concentration in the sample.

# Examples

```julia
myrun = load("data/Lu-Hf",instrument="Agilent")
method = "Lu-Hf"
channels = Dict("d"=>"Hf178 -> 260",
                "D"=>"Hf176 -> 258",
                "P"=>"Lu175 -> 175")
standards = Dict("Hogsbo" => "hogsbo")
glass = Dict("NIST612" => "NIST612p")
blk, fit = process!(myrun,method,channels,standards,glass)
```
"""
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

"""
fitBlanks(run::Vector{Sample};nblank=2):
Fit a dataframe of blank parameters to a run of multiple samples
"""
function fitBlanks(run::Vector{Sample};nblank=2)
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

"""
fractionation

Fit the drift and down hole fractionation

# Methods

- `fractionation(run::Vector{Sample},
                 method::AbstractString,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 standards::Union{AbstractVector,AbstractDict},
                 glass::Union{AbstractVector,AbstractDict};
                 ndrift::Integer=1,
                 ndown::Integer=0,
                 PAcutoff=nothing,
                 verbose::Bool=false)`
- `fractionation(run::Vector{Sample},
                 method::AbstractString,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 standards::Union{AbstractVector,AbstractDict},
                 mf::Union{AbstractFloat,Nothing};
                 ndrift::Integer=1,
                 ndown::Integer=0,
                 PAcutoff=nothing,
                 verbose::Bool=false)`
- `fractionation(run::Vector{Sample},
                 method::AbstractString,
                 blank::AbstractDataFrame,
                 channels::AbstractDict,
                 glass::Union{AbstractVector,AbstractDict};
                 verbose::Bool=false)`
- `fractionation(run::Vector{Sample},
                 blank::AbstractDataFrame,
                 internal::Tuple,
                 glass::Union{AbstractVector,AbstractDict})`
- `fractionation(run::Vector{Sample},
                 blank::AbstractDataFrame,
                 elements::AbstractDataFrame,
                 internal::Tuple,
                 glass::Union{AbstractVector,AbstractDict})`

# Arguments

- see [`process!`](@ref).
- `elements`: a 1-row dataframe with the elements corresponding to each channel
"""
# two-step isotope fractionation
function fractionation(run::Vector{Sample},
                       method::AbstractString,
                       blank::AbstractDataFrame,
                       channels::AbstractDict,
                       standards::AbstractVector,
                       glass::AbstractVector;
                       ndrift::Integer=1,
                       ndown::Integer=0,
                       PAcutoff=nothing,
                       verbose::Bool=false)
    mf = fractionation(run,method,blank,channels,glass;verbose=verbose)
    out = fractionation(run,method,blank,channels,standards,mf;
                        ndrift=ndrift,ndown=ndown,PAcutoff=PAcutoff,verbose=verbose)
    return out
end
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
    return fractionation(run,method,blank,channels,
                         collect(keys(standards)),
                         collect(keys(glass));
                         ndrift=ndrift,ndown=ndown,
                         PAcutoff=PAcutoff,verbose=verbose)
end
# one-step isotope fractionation using mineral standards
function fractionation(run::Vector{Sample},
                       method::AbstractString,
                       blank::AbstractDataFrame,
                       channels::AbstractDict,
                       standards::AbstractVector,
                       mf::Union{AbstractFloat,Nothing};
                       ndrift::Integer=1,
                       ndown::Integer=0,
                       PAcutoff=nothing,
                       verbose::Bool=false)

    anchors = getAnchors(method,standards,false)
    
    if ndrift<1 PTerror("ndriftzero") end

    dats = Dict()
    for (refmat,anchor) in anchors
        dats[refmat] = pool(run;signal=true,group=refmat)
    end

    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]
    bP = blank[:,channels["P"]]

    function misfit(par)
        drift = par[1:ndrift]
        down = vcat(0.0,par[ndrift+1:ndrift+ndown])
        mfrac = isnothing(mf) ? par[ndrift+ndown+1] : log(mf)
        pa = isnothing(PAcutoff) ? 0.0 : par[end]
        out = 0.0
        for (refmat,dat) in dats
            t = dat.t
            T = dat.T 
            Pm = dat[:,channels["P"]]
            Dm = dat[:,channels["D"]]
            dm = dat[:,channels["d"]]
            (x0,y0,y1) = anchors[refmat]
            out += SS(t,T,Pm,Dm,dm,x0,y0,y1,drift,down,mfrac,bP,bD,bd;
                      PAcutoff=PAcutoff,pa=pa)
        end
        return out
    end

    init = fill(0.0,ndrift)
    if (ndown>0) init = vcat(init,fill(0.0,ndown)) end
    if isnothing(mf) init = vcat(init,0.0) end
    if !isnothing(PAcutoff) init = vcat(init,0.0) end
    if length(init)>0
        fit = Optim.optimize(misfit,init)
        if verbose
            println("Drift and downhole fractionation correction:\n")
            println(fit)
        else
            if fit.iteration_converged
                @warn "Reached the maximum number of iterations before achieving " *
                    "convergence. Reduce the order of the polynomials or fix the " *
                    "mass fractionation and try again."
            end
        end
        pars = Optim.minimizer(fit)
    else
        pars = 0.0
    end
    drift = pars[1:ndrift]
    down = vcat(0.0,pars[ndrift+1:ndrift+ndown])
    mfrac = isnothing(mf) ? pars[ndrift+ndown+1] : log(mf)
    pa = isnothing(PAcutoff) ? 0.0 : pars[end]

    return (drift=drift,down=down,mfrac=mfrac,
            PAcutoff=PAcutoff,pa=pa)
end
function fractionation(run::Vector{Sample},
                       method::AbstractString,
                       blank::AbstractDataFrame,
                       channels::AbstractDict,
                       standards::AbstractDict,
                       mf::Union{AbstractFloat,Nothing};
                       ndrift::Integer=1,
                       ndown::Integer=0,
                       PAcutoff=nothing,
                       verbose::Bool=false)
    return fractionation(run,method,blank,channels,
                         collect(keys(standards)),mf;
                         ndrift=ndrift,ndown=ndown,
                         PAcutoff=PAcutoff,verbose=verbose)
end
# isotopic mass fractionation using glass
function fractionation(run::Vector{Sample},
                       method::AbstractString,
                       blank::AbstractDataFrame,
                       channels::AbstractDict,
                       glass::AbstractVector;
                       verbose::Bool=false)
    
    anchors = getAnchors(method,glass,true)

    dats = Dict()
    for (refmat,anchor) in anchors
        dats[refmat] = pool(run;signal=true,group=refmat)
    end

    bD = blank[:,channels["D"]]
    bd = blank[:,channels["d"]]

    function misfit(par)
        mfrac = par[1]
        out = 0.0
        for (refmat,dat) in dats
            t = dat.t
            Dm = dat[:,channels["D"]]
            dm = dat[:,channels["d"]]
            y0 = anchors[refmat]
            out += SS(t,Dm,dm,y0,mfrac,bD,bd)
        end
        return out
    end

    fit = Optim.optimize(misfit,[0.0])
    if verbose
        println("Mass fractionation correction:\n")
        println(fit)
    end

    mfrac = Optim.minimizer(fit)[1]
    
    return exp(mfrac)
end
function fractionation(run::Vector{Sample},
                       method::AbstractString,
                       blank::AbstractDataFrame,
                       channels::AbstractDict,
                       glass::AbstractDict;
                       verbose::Bool=false)
    return fractionation(run,method,blank,channels,
                         collect(keys(glass));
                         verbose=verbose)
end
# for concentration measurements:
function fractionation(run::Vector{Sample},
                       blank::AbstractDataFrame,
                       internal::Tuple,
                       glass::AbstractDict)
    elements = channels2elements(run)
    return fractionation(run,blank,elements,internal,
                         collect(keys(glass)))
end
function fractionation(run::Vector{Sample},
                       blank::AbstractDataFrame,
                       internal::Tuple,
                       glass::AbstractVector)
    elements = channels2elements(run)
    return fractionation(run,blank,elements,internal,glass)
end
function fractionation(run::Vector{Sample},
                       blank::AbstractDataFrame,
                       elements::AbstractDataFrame,
                       internal::Tuple,
                       glass::AbstractVector)
    ne = size(elements,2)
    num = den = fill(0.0,ne-1)
    for SRM in glass
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

"""
atomic(samp::Sample,
       channels::AbstractDict,
       blank::AbstractDataFrame,
       pars::NamedTuple)

# Returns

- `P`, `D`, `d`: Vectors with the inferred 'atomic' parent, daughter and sister signals

# Arguments

See [`process!`](@ref).
"""
function atomic(samp::Sample,
                channels::AbstractDict,
                blank::AbstractDataFrame,
                pars::NamedTuple)
    dat = windowData(samp,signal=true)
    Pm = dat[:,channels["P"]]
    Dm = dat[:,channels["D"]]
    dm = dat[:,channels["d"]]
    ft = polyFac(pars.drift,dat.t)
    FT = polyFac(pars.down,dat.T)
    mf = exp(pars.mfrac)
    bPt = polyVal(blank[:,channels["P"]],dat.t)
    bDt = polyVal(blank[:,channels["D"]],dat.t)
    bdt = polyVal(blank[:,channels["d"]],dat.t)
    D = @. (Dm-bDt)
    P = @. (Pm-bPt)/(ft*FT)
    d = @. (dm-bdt)/mf
    return P, D, d
end
export atomic

"""
averat

Average the 'atomic' isotopic ratios for a sample

# Methods

- `averat(samp::Sample,
          channels::AbstractDict,
          blank::AbstractDataFrame,
          pars::NamedTuple)`
- `averat(run::Vector{Sample},
          channels::AbstractDict,
          blank::AbstractDataFrame,
          pars::NamedTuple;
          PAcutoff=nothing)`

# Returns

- a dataframe of P/D and d/D-ratios with their standard errors and error correlations

# Arguments

See [`process!`](@ref).
"""
function averat(samp::Sample,
                channels::AbstractDict,
                blank::AbstractDataFrame,
                pars::NamedTuple)
    P, D, d = atomic(samp,channels,blank,pars)
    nr = length(P)
    muP = Statistics.mean(P)
    muD = Statistics.mean(D)
    mud = Statistics.mean(d)
    E = Statistics.cov(hcat(P,D,d))
    if mud < 0.0
        mud = 0.0
        E[1,3] = E[3,1] = sum((P.-muP).*d)/(nr-1)
        E[2,3] = E[3,2] = sum((D.-muD).*d)/(nr-1)
        E[3,3] = Statistics.mean(d.^2)
    end
    x = muP/muD
    y = mud/muD
    J = [1/muD -muP/muD^2 0;
         0 -mud/muD^2 1/muD]
    covmat = J * (E/nr) * transpose(J)
    sx = sqrt(covmat[1,1])
    sy = sqrt(covmat[2,2])
    rxy = covmat[1,2]/(sx*sy)
    return [x sx y sy rxy]
end
function averat(run::Vector{Sample},
                channels::AbstractDict,
                blank::AbstractDataFrame,
                pars::NamedTuple;
                PAcutoff=nothing)
    ns = length(run)
    nul = fill(0.0,ns)
    out = DataFrame(name=fill("",ns),x=nul,sx=nul,y=nul,sy=nul,rxy=nul)
    analog = isAnalog(run,channels,PAcutoff)
    for i in 1:ns
        samp = run[i]
        out[i,1] = samp.sname
        if length(pars)==3
            samp_pars = pars
        elseif analog[i]
            samp_pars = pars.analog
        else
            samp_pars = pars.pulse
        end
        out[i,2:end] = averat(samp,channels,blank,samp_pars)
    end
    return out
end
export averat

"""
concentrations

Tabulate chemical concentration data

# Methods

- `concentrations(samp::Sample,
                  blank::AbstractDataFrame,
                  pars::AbstractVector,
                  internal::Tuple)`
- `concentrations(samp::Sample,
                  elements::AbstractDataFrame,
                  blank::AbstractDataFrame,
                  pars::AbstractVector,
                  internal::Tuple)`
- `concentrations(run::Vector{Sample},
                  blank::AbstractDataFrame,
                  pars::AbstractVector,
                  internal::Tuple)`
- `concentrations(run::Vector{Sample},
                  elements::AbstractDataFrame,
                  blank::AbstractDataFrame,
                  pars::AbstractVector,
                  internal::Tuple)`

# Returns

- a dataframe with concentration estimates (in ppm) and their standard errors

# Arguments

- See [`process!`](@ref).
- `elements`: a 1-row dataframe with the elements corresponding to each channel

# Examples
```julia
method = "concentrations"
myrun = load("data/Lu-Hf",instrument="Agilent")
internal = ("Al27 -> 27",1.2e5)
glass = Dict("NIST612" => "NIST612p")
setGroup!(myrun,glass)
blk, fit = process!(myrun,internal,glass;nblank=2)
conc = concentrations(myrun,blk,fit,internal)
```
"""
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
    dat = windowData(samp,signal=true)
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
