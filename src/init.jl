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

function init_methods(csv::AbstractString=joinpath(@__DIR__,"../settings/methods.csv"))
    df = CSV.read(csv, DataFrame)
    out = OrderedDict()
    for row in eachrow(df)
        key = row.method
        value = (P=row.P,D=row.D,d=row.d)
        add2od!(out,key,value)
    end
    return out
end

function init_lambdas(csv::AbstractString=joinpath(@__DIR__,"../settings/lambda.csv"))
    tab = CSV.read(csv, DataFrame)
    out = Dict()
    for row in eachrow(tab)
        out[row.method] = (row["lambda"],row["err"])
    end
    return out
end

function init_iratios(csv::AbstractString=joinpath(@__DIR__,"../settings/iratio.csv"))
    tab = CSV.read(csv, DataFrame)
    out = Dict()
    for row in eachrow(tab)
        isotope = row.isotope
        abundance = row.abundance
        element = row.element
        entry = NamedTuple{(Symbol(isotope),)}((abundance))
        if !(element in keys(out))
            out[element] = entry
        end
        out[element] = merge(out[element],entry)
    end
    return out
end

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

function init_stoichiometry(csv::AbstractString=joinpath(@__DIR__,"../settings/stoichiometry.csv"))
    tab = CSV.read(csv, DataFrame)
    out = OrderedDict()
    good = .!ismissing.(tab)
    for i in axes(tab,1)
        mineral = tab[i,"mineral"]
        conc = tab[i,2:end]
        concvec = collect(conc)
        good = .!ismissing.(concvec)
        add2od!(out,mineral,conc[good])
    end
    return out
end

function init_stoichiometry!(csv::AbstractString=joinpath(@__DIR__,"../settings/stoichiometry.csv"))
    _KJ[]["stoichiometry"] = init_stoichiometry(csv)
end

function init_glass(csv::AbstractString=joinpath(@__DIR__,"../settings/glass.csv"))
    tab = CSV.read(csv, DataFrame)
    out = OrderedDict()
    for row in eachrow(tab)
        add2od!(out,row["SRM"],row[2:end])
    end
    return out
end

function init_glass!(csv::AbstractString=joinpath(@__DIR__,"../settings/glass.csv"))
    _KJ[]["glass"] = init_glass(csv)
end

function init_referenceMaterials(;isochrons::AbstractString=joinpath(@__DIR__,"../settings/standards/isochron.csv"),
                                  points::AbstractString=joinpath(@__DIR__,"../settings/standards/point.csv"),
                                  bias::AbstractString=joinpath(@__DIR__,"../settings/standards/bias.csv"))
    out = Dict()
    for (fname,fun) in Dict(isochrons => IsochronRefmat,
                            points => PointRefmat,
                            bias => BiasRefmat)
        for row in eachrow(CSV.read(fname, DataFrame))
            method = row["method"]
            if !(method in keys(out))
                out[method] = OrderedDict()
            end
            val = fun(row[3],collect(row[4:end-1])...)
            add2od!(out[method],row["name"],val)
        end
    end
    return out
end

function init_referenceMaterials!(;isochrons::AbstractString=joinpath(@__DIR__,"../settings/standards/isochron.csv"),
                                   points::AbstractString=joinpath(@__DIR__,"../settings/standards/point.csv"),
                                   bias::AbstractString=joinpath(@__DIR__,"../settings/standards/bias.csv"))
    _KJ[]["refmat"] = init_referenceMaterials(;isochrons=isochrons,
                                             points=points,
                                             bias=bias)
end
