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

"""
KJmethod(name::String
         channels::DataFrame
         nblanks::Int
         ndrift::Int
         ndown::Int
         standards::Dict
         anchors::Dict)

KJ method type
"""
mutable struct KJmethod
    name::String
    channels::DataFrame
    nblank::Int
    ndrift::Int
    ndown::Int
    standards::Dict
    anchors::Dict
end
export KJmethod

"""
KJfit(blank::DataFrame
      drift::Vector{Real}
      down::Vector{Real}
      PAcutoff::Union{Nothing,Real}
      adrift::Vector{Real})

KJ fit type
"""
mutable struct KJfit
    blank::DataFrame
    drift::Vector{Real}
    down::Vector{Real}
    PAcutoff::Union{Nothing,Real}
    adrift::Vector{Real}
end
export KJfit
