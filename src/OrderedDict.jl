OrderedDict() = OrderedDict(Vector{String}(),Dict())

function add2od!(od::OrderedDict,
                 key::AbstractString,
                 value::Any)
    push!(od.names,String(key))
    push!(od.dict,String(key) => value)
end

function get(od::OrderedDict,
             i::Integer)
    key = od.names[i]
    return get(dict,key)
end
function get(od::OrderedDict,
             key::AbstractString)
    return od.dict[key]
end
