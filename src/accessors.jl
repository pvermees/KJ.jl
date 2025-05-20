"""
getKJctrl

Access the control parameters of a TUI session

# Returns

A dictionary with control parameters

# Examples

See [`KJ`](@ref).

See also [`setKJctrl!`](@ref).
"""
function getKJctrl()
    return _KJ["ctrl"]
end
export getKJctrl

"""
setKJctrl

Set the control parameters of a TUI session

# Examples
```julia
TUI(logbook="logs/test.log")
ctrl = getKJctrl()
ctrl["transformation"] = "log"
setKJctrl!(ctrl)
TUI()
```

See also [`getKJctrl`](@ref).
"""
function setKJctrl!(ctrl::AbstractDict)
    _KJ["ctrl"] = ctrl
end
export setKJctrl!

function getExt(format)
    if format in ["Agilent","ThermoFisher"]
        return ".csv"
    else
        return format
    end
end

function getChannels(run::Vector{Sample})
    return getChannels(run[1])
end
function getChannels(samp::Sample)
    return names(getSignals(samp))
end
export getChannels

function getSnames(run::Vector{Sample})
    return getAttr(run,:sname)
end
export getSnames

function getGroups(run::Vector{Sample})
    return getAttr(run,:group)
end
export getGroups

function getAttr(run::Vector{Sample},
                 attr::Symbol)
    ns = length(run)
    first = getproperty(run[1],attr)
    out = fill(first,ns)
    for i in eachindex(run)
        out[i] = getproperty(run[i],attr)
    end
    return out
end

function setGroup!(run::Vector{Sample},
                   selection::Vector{Int},
                   refmat::AbstractString)
    for i in selection
        run[i].group = refmat
    end
end
function setGroup!(run::Vector{Sample},
                   prefix::AbstractString,
                   refmat::AbstractString)
    snames = getSnames(run)
    selection = findall(contains(prefix),snames)
    setGroup!(run,selection,refmat)
end
function setGroup!(run::Vector{Sample},
                   standards::AbstractDict)
    for (refmat,prefix) in standards
        setGroup!(run,prefix,refmat)
    end
end
function setGroup!(run::Vector{Sample},refmat::AbstractString)
    for sample in run
        sample.group = refmat
    end
end
export setGroup!

function setBwin!(run::Vector{Sample},bwin::AbstractVector;seconds::Bool=false)
    for i in eachindex(run)
        setBwin!(run[i],bwin;seconds=seconds)
    end
end
function setBwin!(samp::Sample,bwin::AbstractVector;seconds::Bool=false)
    samp.bwin = seconds ? time2window(samp,bwin) : bwin
end
function setBwin!(run::Vector{Sample})
    for i in eachindex(run)
        setBwin!(run[i])
    end
end
function setBwin!(samp::Sample)
    bwin = autoWindow(samp,blank=true)
    setBwin!(samp,bwin)
end
export setBwin!

function setSwin!(run::Vector{Sample},swin::AbstractVector;seconds::Bool=false)
    for i in eachindex(run)
        setSwin!(run[i],swin;seconds=seconds)
    end
end
function setSwin!(samp::Sample,swin::AbstractVector;seconds::Bool=false)
    samp.swin = seconds ? time2window(samp,swin) : swin
end
function setSwin!(run::Vector{Sample})
    for samp in run
        setSwin!(samp)
    end
end
function setSwin!(samp::Sample)
    swin = autoWindow(samp,blank=false)
    setSwin!(samp,swin)
end
export setSwin!

function shift_windows!(run::Vector{Sample},
                        shift::Number=0.0)
    for samp in run
        samp.t0 += shift
        di = t2i(samp,shift)
        nt = size(samp.dat,1)
        bwin = []
        for win in samp.bwin
            start = max(win[1]+di,1)
            stop = min(win[2]+di,nt)
            push!(bwin,(start,stop))
        end
        setBwin!(samp,bwin)
        swin = []
        for win in samp.swin
            start = max(win[1]+di,1)
            stop = min(win[2]+di,nt)
            push!(swin,(start,stop))
        end
        setSwin!(samp,swin)
    end
end
export shift_windows!

function geti0(signals::AbstractDataFrame)
    total = sum.(eachrow(signals))
    q = Statistics.quantile(total,[0.05,0.95])
    mid = (q[2]+q[1])/2
    (lovals,lens) = rle(total.<mid)
    i = findfirst(lovals)
    return max(2,sum(lens[1:i]))
end
export geti0

function sett0!(run::Vector{Sample})
    for samp in run
        sett0!(samp)
    end
end
function sett0!(run::Vector{Sample},t0::AbstractFloat)
    for samp in run
        sett0!(samp,t0)
    end
end
function sett0!(samp::Sample)
    dat = getSignals(samp)
    i0 = geti0(dat)
    sett0!(samp,samp.dat[i0,1])
end
function sett0!(samp::Sample,t0::AbstractFloat)
    samp.t0 = t0
end
export sett0!

# isochron
function getx0y0y1(method::AbstractString,
                   refmat::AbstractString)
    t = get(_KJ["refmat"][method],refmat).tx[1]
    if method=="U-Pb"
        L8 = _KJ["lambda"]["U238-Pb206"][1]
        L5 = _KJ["lambda"]["U235-Pb207"][1]
        U58 = _KJ["iratio"]["U-Pb"].U235/_KJ["iratio"]["U-Pb"].U238
        x0 = 1/(exp(L8*t)-1)
        y1 = U58*(exp(L5*t)-1)/(exp(L8*t)-1)
    else
        L = _KJ["lambda"][method][1]
        x0 = 1/(exp(L*t)-1)
        y1 = 0.0
    end
    y0 = get(_KJ["refmat"][method],refmat).y0[1]
    return (x0=x0,y0=y0,y1=y1)
end
# point
function getx0y0(method::AbstractString,
                 refmat::AbstractString)
    x0 = get(_KJ["refmat"][method],refmat).tx[1]
    y0 = get(_KJ["refmat"][method],refmat).y0[1]
    return (x0=x0,y0=y0)
end
# glass
function gety0(method::AbstractString,
               refmat::AbstractString)
    P, D, d = getPDd(method)
    ratio = d * D
    return get(_KJ["glass"],refmat)[ratio]
end

function getAnchors(method::AbstractString,
                    standards::AbstractVector,
                    glass::AbstractVector)
    Ganchors = getGlassAnchors(method,glass)
    Sanchors = getStandardAnchors(method,standards)
    return merge(Ganchors,Sanchors)
end
function getAnchors(method::AbstractString,
                    standards::AbstractDict,
                    glass::AbstractDict)
    return getAnchors(method,collect(keys(standards)),collect(keys(glass)))
end
export getAnchors

function isochronAnchor(anchor::NamedTuple)
    return all(in(keys(anchor)), [:x0,:y0,:y1])
end
function pointAnchor(anchor::NamedTuple)
    k = keys(anchor)
    return in(:x0,k) & in(:y0,k) & !in(:y1,k)
end

function getStandardAnchors(method::AbstractString,
                           refmats::AbstractVector)
    out = Dict()
    for refmat in refmats
        t = get(_KJ["refmat"][method],refmat).type
        if t == "isochron"
            out[refmat] = getx0y0y1(method,refmat)
        else # point
            out[refmat] = getx0y0(method,refmat)
        end
    end
    return out
end
function getStandardAnchors(method::AbstractString,
                           refmats::AbstractDict)
    return getStandardAnchors(method,collect(keys(refmats)))
end
export getStandardAnchors

function getGlassAnchors(method::AbstractString,
                         refmats::AbstractVector,
                         glass::Bool=false)
    out = Dict()
    for refmat in refmats
        out[refmat] = gety0(method,refmat)
    end
    return out
end
function getGlassAnchors(method::AbstractString,
                                refmats::AbstractDict)
    return getGlassAnchors(method,collect(keys(refmats)))
end
export getGlassAnchors

function getSignals(dat::AbstractDataFrame)
    tail = count(x -> x in ["T","x","y"], names(dat)) + 1
    return dat[:,2:end-tail]
end
function getSignals(samp::Sample)
    return getSignals(samp.dat)
end
function getSignals(samp::Sample,channels::AbstractDict)
    return samp.dat[:,collect(values(channels))]
end
export getSignals

function getPDd(method::AbstractString)
    PDd = get(_KJ["methods"],method)
    return PDd.P, PDd.D, PDd.d
end
export getPDd

function getInternal(mineral::AbstractString,channel::AbstractString)
    element = channel2element(channel,collect(keys(_KJ["nuclides"])))
    concentration = get(_KJ["stoichiometry"],mineral)[element[1]] * 1e5
    return (channel,concentration)
end
export getInternal
