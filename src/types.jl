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

_KJ::AbstractDict = Dict()

function init_KJ!()
    _KJ["methods"] = getMethods()
    _KJ["lambda"] = getLambdas()
    _KJ["iratio"] = getiratios()
    _KJ["nuclides"] = getNuclides()
    _KJ["refmat"] = getReferenceMaterials()
    _KJ["glass"] = getGlass()
    _KJ["stoichiometry"] = getStoichiometry()
    _KJ["tree"] = getKJtree()
    _KJ["ctrl"] = nothing
    _KJ["extensions"] = nothing
end
export init_KJ!
