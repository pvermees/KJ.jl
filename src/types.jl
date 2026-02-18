"""
    Sample

Represents a single LA-ICP-MS measurement.

# Fields
- `sname`: Sample name
- `datetime`: Acquisition date and time
- `dat`: DataFrame containing time-series measurements
- `t0`: Time of ablation onset (seconds)
- `bwin`: Vector of blank window tuples (start_index, end_index)
- `swin`: Vector of signal window tuples (start_index, end_index)
- `group`: Group label (e.g., "sample", standard name)
"""
mutable struct Sample
    sname::String
    datetime::DateTime
    dat::DataFrame
    t0::Float64
    bwin::Vector{Tuple}
    swin::Vector{Tuple}
    group::String
end
export Sample

mutable struct OrderedDict
    names::Vector{String}
    dict::Dict
end

"""
    AbstractRefmat

Abstract base type for reference material definitions.

Subtypes include `IsochronRefmat`, `PointRefmat`, and `BiasRefmat`.
"""
abstract type AbstractRefmat end
export AbstractRefmat

mutable struct IsochronRefmat <: AbstractRefmat
    material::String
    t::Float64
    st::Float64
    y0::Float64
    sy0::Float64
end

mutable struct PointRefmat <: AbstractRefmat
    material::String
    x::Float64
    sx::Float64
    y::Float64
    sy::Float64
end

mutable struct BiasRefmat <: AbstractRefmat
    material::String
    y0::Float64
    sy0::Float64
end

"""
    Calibration

Calibration settings for mass bias correction.

# Fields
- `num`: Named tuple with numerator isotope (ion, channel)
- `den`: Named tuple with denominator isotope (ion, channel)
- `standards`: Set of reference material names used for calibration
"""
@kwdef mutable struct Calibration
    num::NamedTuple{(:ion,:channel),Tuple{String,String}} = (ion="",channel="")
    den::NamedTuple{(:ion,:channel),Tuple{String,String}} = (ion="",channel="")
    standards::Set{String} = Set{String}()
end
export Calibration

"""
    AbstractInterference

Abstract base type for interference correction specifications.

Subtypes include `Interference` for general isobaric interferences and
`REEInterference` for rare earth element oxide interferences.
"""
abstract type AbstractInterference end
export AbstractInterference

"""
    Interference <: AbstractInterference

Interference correction specification for isobaric interferences.

# Fields
- `proxy`: Proxy isotope used to estimate the interference
- `channel`: Channel being corrected for interference
- `bias`: Calibration object for any associated bias correction
"""
@kwdef mutable struct Interference <: AbstractInterference
    proxy::String = ""
    channel::String = ""
    bias::Calibration = Calibration()
end
export Interference

"""
    REEInterference <: AbstractInterference

Rare Earth Element (REE) oxide interference correction specification.

# Fields
- `REE`: REE element channel
- `REEO`: REE oxide channel
- `standards`: Set of reference materials used for calibration
"""
@kwdef mutable struct REEInterference <: AbstractInterference
    REE::String = ""
    REEO::String = ""
    standards::Set{String} = Set{String}()
end
export REEInterference

"""
    Pairing

Isotope pairing specification linking ions, proxies, and channels.

# Fields
- `ion`: Isotope notation (e.g., "Pb206")
- `proxy`: Proxy isotope used for calculations (defaults to ion)
- `channel`: Measured channel name (defaults to proxy)
- `interferences`: Dictionary of interference corrections to apply
"""
@kwdef mutable struct Pairing
    ion::String = ""
    proxy::String = ion
    channel::String = proxy
    interferences::Dict{String,AbstractInterference} = Dict{String,Interference}()
end
export Pairing

"""
    AbstractAnchor

Abstract base type for anchor point definitions in isotope ratio space.

Anchor points represent known isotopic compositions of reference materials
used for fractionation correction. Subtypes include `IsochronAnchor`,
`PointAnchor`, and `BiasAnchor`.
"""
abstract type AbstractAnchor end
export AbstractAnchor

mutable struct IsochronAnchor <: AbstractAnchor
    x0::Float64
    y0::Float64
    y1::Float64
end

mutable struct PointAnchor <: AbstractAnchor
    x::Float64
    y::Float64
end

mutable struct BiasAnchor <: AbstractAnchor
    y::Float64
end

"""
    KJmethod

Abstract base type for LA-ICP-MS data reduction methods.

Subtypes include:
- `Gmethod`: Geochronology methods for parent-daughter isotope dating
- `Cmethod`: Concentration methods for quantitative element analysis
"""
abstract type KJmethod end
export KJmethod

"""
    Gmethod <: KJmethod

Geochronology method definition for parent-daughter isotope dating systems.

# Fields
- `name`: Method name (e.g., "U-Pb", "Rb-Sr", "Lu-Hf", "K-Ca", "Re-Os")
- `groups`: Dictionary mapping group names to reference material names
- `P`, `D`, `d`: Pairing definitions for parent, daughter, and normalizing isotopes
- `bias`: Calibration settings for mass bias correction
- `standards`: Set of standard names used for fractionation correction
- `nblank`: Polynomial order for blank fitting (default: 2)
- `ndrift`: Polynomial order for drift correction (default: 2)
- `ndown`: Polynomial order for downhole fractionation (default: 1)
- `nbias`: Polynomial order for bias correction (default: 1)
- `PAcutoff`: Cutoff for analog vs counting mode (default: Inf)
"""
@kwdef mutable struct Gmethod <: KJmethod
    name::String = "U-Pb"
    groups::Dict{String,String} = Dict{String,String}()
    P::Pairing = Pairing(ion=default_ions(name).P)
    D::Pairing = Pairing(ion=default_ions(name).D)
    d::Pairing = Pairing(ion=default_ions(name).d)
    bias::Calibration = Calibration()
    standards::Set{String} = Set(collect(keys(groups)))
    nblank::Int = 2
    ndrift::Int = 2
    ndown::Int = 1
    nbias::Int = 1
    PAcutoff::Float64 = Inf
end
export Gmethod

"""
    Cmethod <: KJmethod

Concentration method definition for quantitative element analysis.

# Fields
- `elements`: Named tuple of elements to quantify
- `groups`: Dictionary mapping group names to reference material names
- `internal`: Tuple specifying internal standard (channel, concentration)
- `nblank`: Polynomial order for blank fitting
"""
mutable struct Cmethod <: KJmethod
    elements::NamedTuple
    groups::Dict{String,String}
    internal::Tuple
    nblank::Int
end
export Cmethod

"""
    AbstractBias

Abstract base type for mass bias correction parameters.

Subtypes include `Bias` for standard isotope ratio bias and `REEBias`
for rare earth element oxide bias corrections.
"""
abstract type AbstractBias end
export AbstractBias

"""
    Bias <: AbstractBias

Mass bias correction parameters for isotope ratios.

# Fields
- `mass_num`: Atomic mass of numerator isotope
- `mass_den`: Atomic mass of denominator isotope
- `par`: Vector of polynomial coefficients for bias correction
"""
mutable struct Bias <: AbstractBias
    mass_num::Int
    mass_den::Int
    par::Vector{Float64}
end

"""
    REEBias <: AbstractBias

Rare Earth Element oxide bias correction parameters.

# Fields
- `par`: Vector of polynomial coefficients for bias correction
"""
mutable struct REEBias <: AbstractBias
    par::Vector{Float64}
end

"""
    KJfit

Abstract base type for fitted LA-ICP-MS correction parameters.

Subtypes include:
- `Gfit`: Fitted parameters for geochronology methods
- `Cfit`: Fitted parameters for concentration methods
"""
abstract type KJfit end
export KJfit

"""
    Gfit <: KJfit

Fitted parameters for geochronology methods.

# Fields
- `blank`: DataFrame of polynomial coefficients for blank corrections
- `drift`: Vector of drift correction parameters
- `down`: Vector of downhole fractionation parameters
- `adrift`: Vector of analog mode drift parameters
- `covmat`: Covariance matrix of fitted parameters
- `bias`: Dictionary of mass bias corrections by element
"""
mutable struct Gfit <: KJfit
    blank::DataFrame
    drift::Vector{Float64}
    down::Vector{Float64}
    adrift::Vector{Float64}
    covmat::Matrix
    bias::Dict{String,AbstractBias}
end
export Gfit

"""
    Cfit <: KJfit

Fitted parameters for concentration methods.

# Fields
- `blank`: DataFrame of polynomial coefficients for blank corrections
- `par`: DataFrame of sensitivity factors for each element
"""
mutable struct Cfit <: KJfit
    blank::DataFrame
    par::DataFrame
end
export Cfit