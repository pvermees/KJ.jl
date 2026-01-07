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

abstract type AbstractRefmat end
export AbstractRefmat

mutable struct IsochronRefmat <: AbstractRefmat
    t::Float64
    st::Float64
    y0::Float64
    sy0::Float64
    material::String
end

mutable struct PointRefmat <: AbstractRefmat
    x::Float64
    sx::Float64
    y::Float64
    sy::Float64
    material::String
end

mutable struct Fractionation
    ions::NamedTuple{(:P,:D,:d)}
    proxies::NamedTuple{(:P,:D,:d)}
    channels::NamedTuple{(:P,:D,:d)}
    standards::Set{String}
    bias::Dict{String,Vector{String}}
end
export Fractionation

mutable struct Interference
    ions::Dict{String,Vector{String}}
    proxies::Dict{String,String}
    channels::Dict{String,String}
    bias::Dict{String,Vector{String}}
end
export Interference

abstract type KJmethod end
export KJmethod

mutable struct Gmethod <: KJmethod
    name::String
    fractionation::Fractionation
    interference::Interference
    nblank::Int
    ndrift::Int
    ndown::Int
    nbias::Int
    PAcutoff::Union{Nothing,Float64}    
end
export Gmethod

mutable struct Cmethod <: KJmethod
    elements::NamedTuple
    standards::Set{String}
    internal::Tuple
    nblank::Int
end
export Cmethod

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

abstract type KJfit end
export KJfit

mutable struct Gfit <: KJfit
    blank::DataFrame
    drift::Vector{Float64}
    down::Vector{Float64}
    adrift::Vector{Float64}
    covmat::Matrix
end
export Gfit

mutable struct Cfit <: KJfit
    blank::DataFrame
    par::DataFrame
end
export Cfit

mutable struct Cruncher
    pmb::Vector{Float64}
    Dombi::Vector{Float64}
    bomb::Vector{Float64}
    bpt::Vector{Float64}
    bDot::Vector{Float64}
    bbot::Vector{Float64}
    vp::Float64
    vD::Float64
    vb::Float64
    spD::Float64
    spb::Float64
    sDb::Float64
    bd::Float64
    t::Vector{Float64}
    T::Vector{Float64}
end
export Cruncher

mutable struct Averager
    Phat::Vector{Float64}
    Dhat::Vector{Float64}
    dhat::Vector{Float64}
    vP::Float64
    vD::Float64
    vd::Float64
    sPD::Float64
    sPd::Float64
    sDd::Float64
end