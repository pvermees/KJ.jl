"""
Sample(sname::String,
       datetime::DateTime,
       dat::DataFrame,
       t0::Float64,
       bwin::Vector{Tuple},
       swin::Vector{Tuple},
       group::String)

The fundamental type underlying all KJ data.
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

abstract type KJmethod end
export KJmethod

mutable struct Gmethod <: KJmethod
    name::String
    channels::DataFrame
    interferences::Dict
    nblank::Int
    ndrift::Int
    ndown::Int
    PAcutoff::Union{Nothing,Float64}
    anchors::Dict
end
export Gmethod

mutable struct Cmethod <: KJmethod
    elements::DataFrame
    refmats::Vector{String}
    concentrations::DataFrame
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
