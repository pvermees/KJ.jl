"""
process!

Fits blanks and fractionation effects

# Returns

- `blank`: a dataframe with blank parameters for all specified channels
- `fit`: a tuple or (if method is omitted) vector of fit parameters

# Methods

- `process!(run::Vector{Sample},
            dt::AbstractDict,
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
- `dt`: the dwell times (in seconds) of the mass channels
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
dt = dwelltime(myrun)
method = "Lu-Hf"
channels = Dict("d"=>"Hf178 -> 260",
                "D"=>"Hf176 -> 258",
                "P"=>"Lu175 -> 175")
standards = Dict("Hogsbo" => "hogsbo")
glass = Dict("NIST612" => "NIST612p")
blk, fit = process!(myrun,dt,method,channels,standards,glass)
```
"""
function process!(run::Vector{Sample},
                  method::AbstractString,
                  channels::AbstractDict,
                  standards::AbstractDict,
                  glass::AbstractDict;
                  dt::Union{AbstractDict,Nothing}=nothing,
                  dead::AbstractFloat=0.0,
                  nblank::Integer=2,ndrift::Integer=1,ndown::Integer=1,
                  PAcutoff=nothing,verbose::Bool=false)
    blank = fitBlanks(run;nblank=nblank)
    setGroup!(run,glass)
    setGroup!(run,standards)
    fit = fractionation(run,method,blank,channels,standards,glass;
                        dt=dt,dead=dead,ndrift=ndrift,ndown=ndown,
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

"""
atomic(samp::Sample,
       dt::AbstractDict,
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
                pars::NamedTuple;
                dt::Union{AbstractDict,Nothing}=nothing,
                dead::AbstractFloat=0.0)
    Pm,Dm,dm,vP,vD,vd,ft,FT,mf,bPt,bDt,bdt =
        SSprep(blank[:,channels["P"]],
               blank[:,channels["D"]],
               blank[:,channels["d"]],
               windowData(samp,signal=true),
               channels,
               pars.mfrac,pars.drift,pars.down;
               dt=dt,dead=dead,
               PAcutoff=pars.PAcutoff,
               adrift=pars.adrift)
    Phat = @. (Pm-bPt)*ft*FT
    Dhat = @. (Dm-bDt)*mf
    dhat = @. (dm-bdt)
    return Phat, Dhat, dhat
end
export atomic

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
