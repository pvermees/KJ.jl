export getKJctrl
""" 
getKJctrl()

Access the control parameters of a TUI session
"""
function getKJctrl()
    return _KJ["ctrl"]
end

export setKJctrl!
"""
setKJctrl!(ctrl::AbstractDict)

Set the control parameters of a TUI session
```
"""
function setKJctrl!(ctrl::AbstractDict)
    _KJ["ctrl"] = ctrl
end

function getExt(format)
    if format in ["Agilent","ThermoFisher"]
        return ".csv"
    else
        return format
    end
end

"""
getChannels(run::Vector{Sample})

Get a vector of column headers of the MS data files
"""
function getChannels(run::Vector{Sample})
    return getChannels(run[1])
end
"""
getChannels(samp::Sample)
"""
function getChannels(samp::Sample)
    return names(getSignals(samp))
end
export getChannels

"""
getSnames(run::Vector{Sample})

Get a vector of sample names
"""
function getSnames(run::Vector{Sample})
    return getAttr(run,:sname)
end
export getSnames

"""
getGroups(run::Vector{Sample})

Get the vector of group names (standards, samples, glasses).
"""
function getGroups(run::Vector{Sample})
    return getAttr(run,:group)
end
export getGroups

"""
getIndicesInGroup(run::Vector{Sample},
                  group::AbstractString)

Get the indices of the samples belonging to 'group'
"""
function getIndicesInGroup(run::Vector{Sample},
                           group::AbstractString)
    return findall(getGroups(myrun) .== group)
end
export getIndicesInGroup

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

"""
setGroup!(run::Vector{Sample},
          selection::Vector{Int},
          refmat::AbstractString)

Change the group allocations for a selection of samples in a run
"""
function setGroup!(run::Vector{Sample},
                   selection::Vector{Int},
                   group::AbstractString)
    for i in selection
        run[i].group = group
    end
end

"""
setGroup!(run::Vector{Sample},
          prefix::AbstractString,
          group::AbstractString)

"""
function setGroup!(run::Vector{Sample},
                   prefix::AbstractString,
                   group::AbstractString)
    snames = getSnames(run)
    selection = findall(contains(prefix),snames)
    setGroup!(run,selection,group)
end

"""
setGroup!(run::Vector{Sample},
          standards::AbstractDict)
"""
function setGroup!(run::Vector{Sample},
                   standards::AbstractDict)
    for (group,prefix) in standards
        setGroup!(run,prefix,group)
    end
end

"""
setGroup!(run::Vector{Sample},
          group::AbstractString)

Assign all the samples to the same group
"""
function setGroup!(run::Vector{Sample},
                   group::AbstractString)
    for sample in run
        sample.group = group
    end
end
export setGroup!

"""
geti0(signals::AbstractDataFrame)

Get the index of 'time zero', i.e. the onset of laser ablation.
"""
function geti0(signals::AbstractDataFrame)
    total = sum.(eachrow(signals))
    q = Statistics.quantile(total,[0.05,0.95])
    mid = (q[2]+q[1])/2
    (lovals,lens) = rle(total.<mid)
    i = findfirst(lovals)
    return max(2,sum(lens[1:i]))
end
export geti0

"""
sett0!(run::Vector{Sample})

Automatically set time zero, i.e. the onset of laser ablation.
"""
function sett0!(run::Vector{Sample})
    for samp in run
        sett0!(samp)
    end
end

"""
sett0!(run::Vector{Sample},
       t0::Number)

Manually set time zero for an entire run.
"""
function sett0!(run::Vector{Sample},
                t0::Number)
    for samp in run
        sett0!(samp,t0)
    end
end

"""
sett0!(samp::Sample)

Automatically set time zero for a single sample.
"""
function sett0!(samp::Sample)
    dat = getSignals(samp)
    i0 = geti0(dat)
    sett0!(samp,samp.dat[i0,1])
end

"""
sett0!(samp::Sample,
       t0::Number)

Manually set time zero for a single sample
"""
function sett0!(samp::Sample,
                t0::Number)
    samp.t0 = t0
end
export sett0!

function get_isochron_anchor(method::AbstractString,
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

function is_isochron_anchor(anchor::NamedTuple)
    return all(in(keys(anchor)), [:x0,:y0,:y1])
end

function get_point_anchor(method::AbstractString,
                          refmat::AbstractString)
    x = get(_KJ["refmat"][method],refmat).tx[1]
    y = get(_KJ["refmat"][method],refmat).y0[1]
    return (x=x,y=y)
end

function is_point_anchor(anchor::NamedTuple)
    k = keys(anchor)
    return all(in(keys(anchor)), [:x,:y])
end

function get_glass_anchor(method::AbstractString,
                          refmat::AbstractString)
    P, D, d = getPDd(method)
    ratio = d * D
    y = get(_KJ["glass"],refmat)[ratio]
    return (y=y)
end

"""
getAnchors(method::AbstractString,
           standards::AbstractVector,
           glass::AbstractVector)

Returns a Dict with named tuple containing
x0, y0, y1, x or y for the standards and glasses
"""
function getAnchors(method::AbstractString,
                    standards::AbstractVector,
                    glass::AbstractVector)
    Ganchors = getGlassAnchors(method,glass)
    Sanchors = getStandardAnchors(method,standards)
    return merge(Ganchors,Sanchors)
end

"""
getAnchors(method::AbstractString,
           standards::AbstractDict,
           glass::AbstractDict)
"""
function getAnchors(method::AbstractString,
                    standards::AbstractDict,
                    glass::AbstractDict)
    return getAnchors(method,collect(keys(standards)),collect(keys(glass)))
end
export getAnchors

"""
getStandardAnchors(method::AbstractString,
                   refmats::AbstractVector)

Returns a Dict with named tuple containing
x0, y0, y1, x or y for the age standards
"""
function getStandardAnchors(method::AbstractString,
                            refmats::AbstractVector)
    out = Dict()
    for refmat in refmats
        t = get(_KJ["refmat"][method],refmat).type
        if t == "isochron"
            out[refmat] = get_isochron_anchor(method,refmat)
        else
            out[refmat] = get_point_anchor(method,refmat)
        end
    end
    return out
end

"""
getStandardAnchors(method::AbstractString,
                   refmats::AbstractDict)
"""
function getStandardAnchors(method::AbstractString,
                            refmats::AbstractDict)
    return getStandardAnchors(method,collect(keys(refmats)))
end
export getStandardAnchors

"""
getGlassAnchors(method::AbstractString,
                refmats::AbstractVector)

Returns a Dict with named tuple containing
x0, y0, y1, x or y for the reference glasses
"""
function getGlassAnchors(method::AbstractString,
                         refmats::AbstractVector)
    out = Dict()
    for refmat in refmats
        out[refmat] = get_glass_anchor(method,refmat)
    end
    return out
end

"""
getGlassAnchors(method::AbstractString,
                refmats::AbstractDict)
"""
function getGlassAnchors(method::AbstractString,
                         refmats::AbstractDict)
    return getGlassAnchors(method,collect(keys(refmats)))
end
export getGlassAnchors

function getSignals(dat::AbstractDataFrame)
    tail = count(x -> x in ["outlier","t","T","x","y"], names(dat))
    return dat[:,2:end-tail]
end

"""
getSignals(samp::Sample)

Returns a dataframe with signals (no time, coordinates or outliers)
"""
function getSignals(samp::Sample)
    return getSignals(samp.dat)
end

"""
getSignals(samp::Sample,
           channels::AbstractDict)
"""
function getSignals(samp::Sample,
                    channels::AbstractDict)
    return samp.dat[:,collect(values(channels))]
end
export getSignals

"""
getPDd(method::AbstractString)

Get the names of the parent, daughter and sister channel
"""
function getPDd(method::AbstractString)
    PDd = get(_KJ["methods"],method)
    return PDd.P, PDd.D, PDd.d
end
export getPDd

"""
getInternal(mineral::AbstractString,
            channel::AbstractString)

Get a tuple with the channel and its reference concentration
"""
function getInternal(mineral::AbstractString,
                     channel::AbstractString)
    element = channel2element(channel,collect(keys(_KJ["nuclides"])))
    concentration = get(_KJ["stoichiometry"],mineral)[element[1]] * 1e5
    return (channel,concentration)
end
export getInternal

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
