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

mutable struct Channels
    P::String
    D::String
    d::String
    S::String
end

function Channels(;
                  P::AbstractString="",
                  D::AbstractString="",
                  d::AbstractString="",
                  S::AbstractString="")
    return Channels(P,D,d,S)
end
