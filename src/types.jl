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

@kwdef mutable struct Interference <: AbstractInterference
    ion::String = ""
    proxy::String = ""
    channel::String = ""
    bias_key::String = ""
end
export Interference

@kwdef mutable struct REEInterference <: AbstractInterference
    channel::String = ""
    bias_key::String = ""
end
export REEInterference

@kwdef mutable struct Pairing
    ion::String = ""
    proxy::String = ion
    channel::String = proxy
    interferences::Set = Set{AbstractInterference}()
end
export Pairing

@kwdef mutable struct Calibration
    num::NamedTuple{(:ion,:channel),Tuple{String,String}} = (ion="",channel="")
    den::NamedTuple{(:ion,:channel),Tuple{String,String}} = (ion="",channel="")
    standards::Set{String} = Set{String}()
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

mutable struct BiasAnchor <: AbstractAnchor
    y::Float64
end

abstract type KJmethod end
export KJmethod

@kwdef mutable struct Gmethod <: KJmethod
    name::String = "U-Pb"
    groups::Dict{String,String} = Dict{String,String}()
    P::Pairing = Pairing(ion=default_ions(name).P)
    D::Pairing = Pairing(ion=default_ions(name).D)
    d::Pairing = Pairing(ion=default_ions(name).d)
    bias::Dict = Dict{String,Calibration}()
    standards::Set{String} = Set(collect(keys(groups)))
    nblank::Int = 2
    ndrift::Int = 2
    ndown::Int = 1
    nbias::Int = 1
    PAcutoff::Float64 = Inf
end
export Gmethod

mutable struct Cmethod <: KJmethod
    elements::NamedTuple
    groups::Dict{String,String}
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