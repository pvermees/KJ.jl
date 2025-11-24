"""
PDdS(;
     P::String="",
     D::String="",
     d::String="",
     S::Union{Missing,String}=missing)
"""
function PDdS(;
              P::String="",
              D::String="",
              d::String="",
             S::Union{Missing,String}=missing)
    return PDdS(P,D,d,S)
end
"""
PDdS(method::AbstractString;
     P::Union{String,Nothing}=nothing,
     D::Union{String,Nothing}=nothing,
     d::Union{String,Nothing}=nothing,
     S::Union{String,Nothing}=nothing)
"""
function PDdS(method::AbstractString;
              P::Union{String,Nothing}=nothing,
              D::Union{String,Nothing}=nothing,
              d::Union{String,Nothing}=nothing,
              S::Union{String,Nothing}=nothing)
    PDdS = get(_KJ["methods"],method)
    return PDdS(ifelse(isnothing(P),PDdS.P,P),
                ifelse(isnothing(D),PDdS.D,D),
                ifelse(isnothing(d),PDdS.d,d),
                ifelse(isnothing(S),PDdS.S,S))
end
