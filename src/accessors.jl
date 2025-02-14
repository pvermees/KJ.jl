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

function getExt(instrument)
    if instrument in ["Agilent","ThermoFisher"]
        return ".csv"
    else
        KJerror("unknownInstrument")
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
    for i in eachindex(run)
        setSwin!(run[i])
    end
end
function setSwin!(samp::Sample)
    swin = autoWindow(samp,blank=false)
    setSwin!(samp,swin)
end
export setSwin!

function geti0(signals::AbstractDataFrame)
    total = sum.(eachrow(signals))
    q = Statistics.quantile(total,[0.05,0.95])
    mid = (q[2]+q[1])/2
    (lovals,lens) = rle(total.<mid)
    i = findfirst(lovals)
    return sum(lens[1:i])
end
export geti0

function sett0!(samp::Sample)
    dat = getSignals(samp)
    i0 = geti0(dat)
    samp.t0 = samp.dat[i0,1]
end
export sett0!

# mineral
function getx0y0y1(method::AbstractString,
                   refmat::AbstractString)
    t = _KJ["refmat"][method][refmat].t[1]
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
    y0 = _KJ["refmat"][method][refmat].y0[1]
    return (x0=x0,y0=y0,y1=y1)
end
# glass
function gety0(method::AbstractString,
               refmat::AbstractString)
    i = findfirst(==(method),_KJ["methods"][:,"method"])
    ratio = _KJ["methods"][i,"d"] * _KJ["methods"][i,"D"]
    return _KJ["glass"][refmat][ratio]
end

function getAnchors(method::AbstractString,
                    standards::AbstractVector,
                    glass::AbstractVector)
    Sanchors = getAnchors(method,standards,false)
    Ganchors = getAnchors(method,glass,true)
    return merge(Sanchors,Ganchors)
end
function getAnchors(method::AbstractString,
                    standards::AbstractDict,
                    glass::AbstractDict)
    return getAnchors(method,collect(keys(standards)),collect(keys(glass)))
end
function getAnchors(method::AbstractString,
                    refmats::AbstractVector,
                    glass::Bool=false)
    out = Dict()
    for refmat in refmats
        out[refmat] = glass ? gety0(method,refmat) : getx0y0y1(method,refmat)
    end
    return out
end
function getAnchors(method::AbstractString,
                    refmats::AbstractDict,
                    glass::Bool=false)
    return getAnchors(method,collect(keys(refmats)),glass)
end
export getAnchors

function getSignals(dat::AbstractDataFrame)
    tail = "T" in names(dat) ? 2 : 1
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
    i = findfirst(==(method),_KJ["methods"][:,"method"])
    PDd = _KJ["methods"][i,2:end]
    return PDd.P, PDd.D, PDd.d
end
export getPDd
function getMethods(csv::AbstractString=joinpath(@__DIR__,"../settings/methods.csv"))
    return CSV.read(csv, DataFrame)
end
export getMethods
function getLambdas(csv::AbstractString=joinpath(@__DIR__,"../settings/lambda.csv"))
    tab = CSV.read(csv, DataFrame)
    out = Dict()
    for row in eachrow(tab)
        out[row.method] = (row["lambda"],row["err"])
    end
    return out
end
export getLambdas
function getiratios(csv::AbstractString=joinpath(@__DIR__,"../settings/iratio.csv"))
    tab = CSV.read(csv, DataFrame)
    out = Dict()
    for row in eachrow(tab)
        isotope = row.isotope
        abundance = row.abundance
        method = row.method
        entry = NamedTuple{(Symbol(isotope),)}((abundance))
        if !(method in keys(out))
            out[method] = entry
        end
        out[method] = merge(out[method],entry)
    end
    return out
end
export getiratios
function getNuclides(csv::AbstractString=joinpath(@__DIR__,"../settings/nuclides.csv"))
    tab = CSV.read(csv, DataFrame)
    elements = unique(tab[:,:element])
    out = Dict()
    for element in elements
        i = findall(tab[:,:element] .== element)
        out[element] = tab[i,:isotope]
    end
    return out
end
export getNuclides
function getStoichiometry(csv::AbstractString=joinpath(@__DIR__,"../settings/stoichiometry.csv"))
    tab = CSV.read(csv, DataFrame)
    out = Dict()
    good = .!ismissing.(tab)
    (nr,nc) = size(tab)
    for i in 1:nr
        mineral = tab[i,"mineral"]
        out[mineral] = Dict()
        for j in 2:nc
            element = names(tab)[j]
            concentration = tab[i,j]
            if !ismissing(concentration)
                out[mineral][element] = concentration
            end
        end
    end
    return out
end
export getStoichiometry
function setStoichiometry(csv::AbstractString=joinpath(@__DIR__,"../settings/stoichiometry.csv"))
    _KJ["stoichiometry"] = getStoichiometry(csv)
end
export setStoichiometry
function getGlass(csv::AbstractString=joinpath(@__DIR__,"../settings/glass.csv"))
    tab = CSV.read(csv, DataFrame)
    out = Dict()
    for row in eachrow(tab)
        out[row["SRM"]] = row[2:end]
    end
    return out
end
export getGlass
function setGlass!(csv::AbstractString=joinpath(@__DIR__,"../settings/glass.csv"))
    _KJ["glass"] = getGlass(csv)
end
function getReferenceMaterials(csv::AbstractString=joinpath(@__DIR__,"../settings/standards.csv"))
    tab = CSV.read(csv, DataFrame)
    out = Dict()
    for row in eachrow(tab)
        method = row["method"]
        if !(method in keys(out))
            out[method] = Dict()
        end
        name = row["name"]
        out[method][name] = (t=(row["t"],row["st"]),y0=(row["y0"],row["sy0"]))
    end
    return out
end
export getReferenceMaterials
function setReferenceMaterials!(csv::AbstractString=joinpath(@__DIR__,"../settings/standards.csv"))
    _KJ["refmat"] = getReferenceMaterials(csv)
end
export setReferenceMaterials!
function getInternal(mineral::AbstractString,channel::AbstractString)
    element = channel2element(channel,collect(keys(_KJ["nuclides"])))
    concentration = _KJ["stoichiometry"][mineral][element[1]] * 1e5
    return (channel,concentration)
end
export getInternal
