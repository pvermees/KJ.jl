OrderedDict() = OrderedDict(Vector{String}(),Dict())

function add2od!(od::OrderedDict,
                 key::AbstractString,
                 value::Any)
    push!(od.names,String(key))
    push!(od.dict,String(key) => value)
end

function get(od::OrderedDict,
             i::Integer)
    n = length(od.names)
    key = od.names[mod1(i,n)]
    return get(od,key)
end
function get(od::OrderedDict,
             key::AbstractString)
    return od.dict[string(key)]
end
