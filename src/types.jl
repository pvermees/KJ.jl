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
PDdS(P::String,
     D::String,
     d::String,
     S::Union{Missing,String})

P = Parent, D = Daughter, d = sister, S = Stranger
"""
mutable struct PDdS
    P::String
    D::String
    d::String
    S::Union{Missing,String}
end
export PDdS
