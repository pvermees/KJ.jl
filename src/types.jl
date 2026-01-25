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

abstract type AbstractInterference end
export AbstractInterference

mutable struct Interference <: AbstractInterference
    ion::String
    proxy::String
    channel::String
    standards::Set{String}
end

mutable struct REEInterference <: AbstractInterference
    proxychannel::String
    numchannel::String
    denchannel::String
end

mutable struct Pairing
    ion::String
    proxy::String
    channel::String
    interferences::Set{AbstractInterference}
end
export Pairing

mutable struct Calibration
    num::NamedTuple{(:ion,:channel),Tuple{String,String}}
    den::NamedTuple{(:ion,:channel),Tuple{String,String}}
    standards::Set{String}
end
export Calibration

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

abstract type KJmethod end
export KJmethod

mutable struct Gmethod <: KJmethod
    name::String
    groups::Dict{String,String}
    P::Pairing
    D::Pairing
    d::Pairing
    bias::Dict{String,Calibration}
    standards::Set{String}
    nblank::Int
    ndrift::Int
    ndown::Int
    nbias::Int
    PAcutoff::Float64
end
export Gmethod

mutable struct Cmethod <: KJmethod
    elements::NamedTuple
    standards::Set{String}
    internal::Tuple
    nblank::Int
end
export Cmethod

abstract type KJfit end
export KJfit

mutable struct Gfit <: KJfit
    blank::DataFrame
    drift::Vector{Float64}
    down::Vector{Float64}
    adrift::Vector{Float64}
    covmat::Matrix
    bias::DataFrame
end
export Gfit

mutable struct Cfit <: KJfit
    blank::DataFrame
    par::DataFrame
end
export Cfit