"""
    concentrations(samp::Sample, method::Cmethod, fit::Cfit; internal=method.internal)
    concentrations(run::Vector{Sample}, method::Cmethod, fit::Cfit)

Calculate element concentrations from calibrated LA-ICP-MS data.

For single samples, returns time-resolved concentrations for each element.
For runs, returns summary statistics (mean ± standard error) for each sample.

# Arguments
- `samp`/`run`: Sample or vector of samples
- `method`: Concentration method with element definitions and internal standard
- `fit`: Fitted sensitivity factors
- `internal`: Internal standard as (channel, concentration) tuple

# Returns
- DataFrame with concentration values. For single samples, includes time-resolved
  data and optional x,y coordinates. For runs, includes sample names, means, and
  standard errors.
"""
function concentrations(samp::Sample,
    method::Cmethod,
    fit::Cfit;
    internal::Dict{String,Tuple{String,Real}}=method.internal
)
    dat = swinData(samp; add_xy=true)
    sig = getSignals(dat)
    bt = predict(samp, fit.blank; t=dat.t)
    X = sig .- bt
    isname, Cs = get_internal(samp, internal)
    Xs = X[:, isname]
    out = (X .* Cs) ./ (Xs .* fit.par)
    nms = "ppm[" .* collect(string.(values(method.elements))) .* "] from " .* names(sig)
    if "x" in names(dat) && "y" in names(dat)
        out.x = dat.x
        out.y = dat.y
        append!(nms, ["x", "y"])
    end
    rename!(out, Symbol.(nms))
    return out
end
function concentrations(run::Vector{Sample},
    method::Cmethod,
    fit::Cfit)
    nr = length(run)
    ne = length(method.elements)
    nc = 2*ne
    mat = zeros(nr, nc)
    col_names = nothing
    for (key, refmat) in method.groups
        if !haskey(method.internal, key)
            isname, _ = get_internal("default", method.internal)
            refconcs = getConcentrations(method, refmat)
            push!(method.internal, Pair(method.groups[key], (isname, refconcs[1, isname])))
        end
    end
    Threads.@threads for i in eachindex(run)
        conc = concentrations(run[i], method, fit)[:, 1:ne]
        if i == firstindex(run)
            col_names = names(conc)
        end
        conc_mat = Matrix(conc)

        mu = vec(Statistics.mean(conc_mat, dims=1))
        sigma = vec(Statistics.std(conc_mat, dims=1))

        nt = size(conc_mat, 1)
        mat[i, 1:2:(nc-1)] .= mu
        mat[i, 2:2:nc] .= sigma ./ sqrt(nt)
    end
    nms = fill("", nc)
    nms[1:2:(nc-1)] .= col_names
    nms[2:2:nc] .= "s[" .* col_names .* "]"
    out = hcat(DataFrame(name=getSnames(run)), DataFrame(mat, Symbol.(nms)))
    return out
end
export concentrations


function get_internal(samp, internal)

    nm = typeof(samp) <: AbstractString ? samp : samp.sname

    for (k, v) in internal
        k == "default" && continue
        occursin(k, nm) && return v
    end

    return internal["default"]
end
