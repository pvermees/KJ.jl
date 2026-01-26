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

function getChannels(run::Vector{Sample}) :: AbstractVector
    return getChannels(run[1])
end
function getChannels(samp::Sample) :: AbstractVector
    return names(getSignals(samp))
end
function getChannels(method::Gmethod) :: AbstractVector
    return [method.P.channel;method.D.channel;method.d.channel]
end
function getChannels(method::Cmethod) :: AbstractVector
    return collect(string.(keys(method.elements)))
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

function getIndicesInGroup(run::Vector{Sample},
                           group::AbstractString)
    return findall(getGroups(run) .== group)
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

function setGroup!(run::Vector{Sample};
                   selection=eachindex(run),
                   group::AbstractString="sample")
    for i in selection
        run[i].group = group
    end
end
function setGroup!(run::Vector{Sample},
                   prefixes::Vector{String})
    snames = getSnames(run)
    for prefix in prefixes
        setGroup!(run;
                  selection=findall(s -> contains(s,prefix), snames),
                  group=prefix)
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
    element = channel2element(channel)
    concentration = get(_KJ["stoichiometry"],mineral)[element] * 1e5
    return (channel,concentration)
end
export getInternal

function getStandards(bias::AbstractDict)
    out = String[]
    for standards in values(bias)
        append!(out,standards)
    end
    return out
end