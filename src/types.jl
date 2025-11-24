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
Chronometer(method::String
            channels::DataFrame
            nblanks::Int
            ndrift::Int
            ndown::Int
            PAcutoff::Union{Nothing,Real})

ions, proxies, channels and interferences are
NamedTuples with keys P, D d and S
"""
mutable struct Chronometer
    method::String
    channels::DataFrame
    nblank::Int
    ndrift::Int
    ndown::Int
    PAcutoff::Union{Nothing,Real}
end
export Chronometer

"""
Fit(blank::DataFrame
    drift::Vector{Real}
    down::Vector{Real}
    adrift::Vector{Real})

KJ standard fit parameters
"""
mutable struct Fit
    blank::DataFrame
    drift::Vector{Real}
    down::Vector{Real}
    adrift::Vector{Real}
end
export Fit
