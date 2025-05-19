function init_KJ!()
    _KJ["methods"] = init_methods()
    _KJ["lambda"] = init_lambdas()
    _KJ["iratio"] = init_iratios()
    _KJ["nuclides"] = init_nuclides()
    _KJ["refmat"] = init_referenceMaterials()
    _KJ["glass"] = init_glass()
    _KJ["stoichiometry"] = init_stoichiometry()
    _KJ["tree"] = init_KJtree()
    _KJ["ctrl"] = nothing
    _KJ["extensions"] = nothing
end
export init_KJ!

function init_methods(csv::AbstractString=joinpath(@__DIR__,"../settings/methods.csv"))
    return CSV.read(csv, DataFrame)
end
export init_methods

function init_lambdas(csv::AbstractString=joinpath(@__DIR__,"../settings/lambda.csv"))
    tab = CSV.read(csv, DataFrame)
    out = Dict()
    for row in eachrow(tab)
        out[row.method] = (row["lambda"],row["err"])
    end
    return out
end
export init_lambdas

function init_iratios(csv::AbstractString=joinpath(@__DIR__,"../settings/iratio.csv"))
    tab = CSV.read(csv, DataFrame)
    out = Dict()
    for row in eachrow(tab)
        isotope = row.isotope
        abundance = row.abundance
        method = row.method
        entry = NamedTuple{(Symbol(isotope),)}((abundance))
        if !(method in keys(out))
            out[method] = entry
        end
        out[method] = merge(out[method],entry)
    end
    return out
end
export init_iratios

function init_nuclides(csv::AbstractString=joinpath(@__DIR__,"../settings/nuclides.csv"))
    tab = CSV.read(csv, DataFrame)
    elements = unique(tab[:,:element])
    out = Dict()
    for element in elements
        i = findall(tab[:,:element] .== element)
        out[element] = tab[i,:isotope]
    end
    return out
end
export init_nuclides

function init_stoichiometry(csv::AbstractString=joinpath(@__DIR__,"../settings/stoichiometry.csv"))
    tab = CSV.read(csv, DataFrame)
    out = Dict()
    good = .!ismissing.(tab)
    (nr,nc) = size(tab)
    for i in 1:nr
        mineral = tab[i,"mineral"]
        out[mineral] = Dict()
        for j in 2:nc
            element = names(tab)[j]
            concentration = tab[i,j]
            if !ismissing(concentration)
                out[mineral][element] = concentration
            end
        end
    end
    return out
end
export init_stoichiometry

function init_stoichiometry!(csv::AbstractString=joinpath(@__DIR__,"../settings/stoichiometry.csv"))
    _KJ["stoichiometry"] = init_stoichiometry(csv)
end
export init_stoichiometry!

function init_glass(csv::AbstractString=joinpath(@__DIR__,"../settings/glass.csv"))
    tab = CSV.read(csv, DataFrame)
    out = Dict()
    for row in eachrow(tab)
        out[row["SRM"]] = row[2:end]
    end
    return out
end
export init_glass

function init_glass!(csv::AbstractString=joinpath(@__DIR__,"../settings/glass.csv"))
    _KJ["glass"] = init_glass(csv)
end
export init_glass!

function init_referenceMaterials(csv::AbstractString=joinpath(@__DIR__,"../settings/standards.csv"))
    tab = CSV.read(csv, DataFrame)
    out = Dict()
    for row in eachrow(tab)
        method = row["method"]
        if !(method in keys(out))
            out[method] = Dict()
        end
        name = row["name"]
        out[method][name] = (tx=(row["tx"],row["stx"]),y0=(row["y0"],row["sy0"]),type=row["type"])
    end
    return out
end
export init_referenceMaterials

function init_referenceMaterials!(csv::AbstractString=joinpath(@__DIR__,"../settings/standards.csv"))
    _KJ["refmat"] = init_referenceMaterials(csv)
end
export init_referenceMaterials!

function init_referenceMaterials!(csv::AbstractString=joinpath(@__DIR__,"../settings/standards.csv"))
    _KJ["refmat"] = init_referenceMaterials(csv)
end
export init_referenceMaterials!
