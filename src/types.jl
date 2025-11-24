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
Method(name::String,
       ions::NamedTuple,
       proxies::Union{Nothing,NamedTuple},
       channels::Union{Nothing,NamedTuple})

ions, proxies, channels and interferences are
NamedTuples with keys P, D d and S
"""
mutable struct Method
    name::String
    ions::NamedTuple,
    proxies::Union{Nothing,NamedTuple}
    channels::Union{Nothing,NamedTuple}
end
export Method
