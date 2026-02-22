abstract type AbstractCache end

mutable struct ProxySelectionCache <: AbstractCache
    isotopes::Vector{String}
    channels::Vector{String}
end

function ProxySelectionCache(;ions::Vector{String}=String[], 
                              channels::Vector{String}=String[],
                              isotopes::Vector{String}=channel2proxy.(ions))
    return ProxySelectionCache(isotopes, channels)
end

@kwdef mutable struct InterferenceCache <: AbstractCache
    target::Pairing = Pairing()
    key::String = ""
    interference::AbstractInterference = Interference()
    glass::String = ""
end

@kwdef mutable struct BiasCache <: AbstractCache
    element::String = ""
    bias::Calibration = Calibration()
    standard::String = ""
end

function add_glass_to_cache(cache::Any,
                            glass::String)
    return glass
end
function add_glass_to_cache(cache::InterferenceCache,
                            glass::String)
    cache.glass = glass
    return cache
end
function add_glass_to_cache(cache::BiasCache,
                            glass::String)
    cache.standard = glass
    return cache
end

function get_glass_from_cache(cache::Any)
    return cache
end
function get_glass_from_cache(cache::InterferenceCache)
    return cache.glass
end
function get_glass_from_cache(cache::BiasCache)
    return cache.standard
end

function get_refmat_from_cache(cache::Any)
    return ""
end
function get_refmat_from_cache(cache::String)
    return cache
end
function get_refmat_from_cache(cache::InterferenceCache)
    return cache.glass
end
function get_refmat_from_cache(cache::BiasCache)
    return cache.standard
end

function push_glass_to_cache!(cache::Any,
                              response::String;
                              ctrl::AbstractDict)
    ctrl["method"].groups[response] = cache
end
function push_glass_to_cache!(cache::InterferenceCache,
                              response::String;
                              other...)
    push!(cache.interference.standards,response)
    interferences = cache.target.interferences
    interferences[cache.key] = cache.interference
end
function push_glass_to_cache!(cache::BiasCache,
                              response::String;
                              ctrl::AbstractDict)
    push!(cache.bias.standards,response)
    TUIaddBias2method!(ctrl)
end

function push_standard_to_cache!(cache::BiasCache,
                                 response::String,
                                 ctrl::AbstractDict)
    standard = cache.standard
    push!(cache.bias.standards,response)
    ctrl["method"].groups[response] = standard
    return TUIaddBias2method!(ctrl)
end
function push_standard_to_cache!(cache::Any,
                                 response::String,
                                 ctrl::AbstractDict)
    if ctrl["method"] isa Gmethod
        push!(ctrl["method"].standards,response)
    end
    ctrl["method"].groups[response] = cache
    setGroup!(ctrl["run"],ctrl["method"])
    ctrl["priority"]["fractionation"] = false
    return "xxx"
end
function push_standard_to_cache!(cache::String,
                                 selection::AbstractVector,
                                 ctrl::AbstractDict)
    setGroup!(ctrl["run"],selection,cache)
    ctrl["method"].groups[cache] = cache
    if ctrl["method"] isa Gmethod
        push!(ctrl["method"].standards,cache)
    end
    ctrl["priority"]["fractionation"] = false
    return "xxx"
end