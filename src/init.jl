_KJ::AbstractDict = Dict()

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
    df = CSV.read(csv, DataFrame)
    out = OrderedDict()
    for row in eachrow(df)
        key = row[1]
        value = (P=String(row[2]),
                 D=String(row[3]),
                 d=String(row[4]))
        add2od!(out,key,value)
    end
    return out
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
    out = OrderedDict()
    good = .!ismissing.(tab)
    (nr,nc) = size(tab)
    for i in 1:nr
        mineral = tab[i,"mineral"]
        conc = tab[i,2:end]
        concvec = collect(conc)
        good = .!ismissing.(concvec)
        add2od!(out,mineral,conc[good])
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
    out = OrderedDict()
    for row in eachrow(tab)
        add2od!(out,row["SRM"],row[2:end])
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
            out[method] = OrderedDict()
        end
        val = (tx=(row["tx"],row["stx"]),
               y0=(row["y0"],row["sy0"]),
               type=row["type"])
        add2od!(out[method],row["name"],val)
    end
    return out
end
export init_referenceMaterials

function init_referenceMaterials!(csv::AbstractString=joinpath(@__DIR__,"../settings/standards.csv"))
    _KJ["refmat"] = init_referenceMaterials(csv)
end
export init_referenceMaterials!
