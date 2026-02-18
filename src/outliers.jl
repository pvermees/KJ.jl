"""
    detect_outliers(vec::AbstractVector; b=2)
    detect_outliers(mat::Matrix; b=2)

Detect outliers in time-series data using a moving median approach.

For vectors, uses interquartile range (IQR) criterion on deviations from moving median.
For matrices, uses Mahalanobis distance and minimum covariance determinant (MCD).

# Arguments
- `vec`/`mat`: Time-series data as vector or matrix
- `b`: Bandwidth for moving median (default: 2)

# Returns
- Vector of indices where outliers are detected
"""
function detect_outliers(vec::AbstractVector;
                         b::Int=2)
    n = length(vec)
    i = moving_median_indices(n;b=b)
    pred = moving_median(vec,i)
    d = pred .- vec
    Q1 = Statistics.quantile(d,0.25)
    Q3 = Statistics.quantile(d,0.75)
    IQR = Q3 - Q1
    ll = Q1 - 1.5 * IQR
    ul = Q3 + 1.5 * IQR
    outliers = @. d > ul || d < ll
    return findall(outliers)
end

function detect_outliers(mat::Matrix;
                         b::Int=2)
    n = size(mat,1)
    i = moving_median_indices(n;b=b)
    pred = moving_median(mat,i)
    d = pred .- mat
    mcd = fast_mcd(d, num_starts=100)
    Sigma = mcd.scatter_raw
    M = d * inv(Sigma) * transpose(d)
    maha = sqrt.(diag(M))
    return detect_outliers(vec(maha))
end
export detect_outliers

"""
    detect_outliers!(run::Vector{Sample}; channels=getChannels(run), include_samples=false)
    detect_outliers!(run::Vector{Sample}, method::KJmethod; include_samples=false)

Detect and flag outliers in a run of samples (in-place modification).

Sets the `outlier` column in each sample's data frame to true for detected outliers.

# Arguments
- `run`: Vector of samples to process
- `channels`/`method`: Channels to use for outlier detection, or method from which to extract channels
- `include_samples`: If true, also detect outliers in samples (not just standards)
"""
function detect_outliers!(run::Vector{Sample};
                          channels::AbstractVector=getChannels(run),
                          include_samples::Bool=false)
    for samp in run
        samp.dat.outlier = falses(size(samp.dat,1))
        if samp.group !== "sample" || include_samples
            blk = bwinData(samp)
            win2outliers!(samp,blk[:,channels],:bwin)
            sig = swinData(samp)
            win2outliers!(samp,sig[:,channels],:swin)
        end
    end
end
function detect_outliers!(run::Vector{Sample}, 
                          method::KJmethod;
                          include_samples::Bool=false)
    ch = getChannels(method)
    detect_outliers!(run;channels=ch,include_samples=include_samples)
end
export detect_outliers!

function win2outliers!(samp::Sample,
                       df::AbstractDataFrame,
                       window::Symbol)
    i = []
    for win in getfield(samp,window)
        append!(i,collect(win[1]:win[2]))
    end
    M = Matrix(df)
    if all(M .> 0)
        L = log.(M)
        x = L[:,1] .- L[:,2:end]
    else
        x = sum(eachcol(M))
    end
    outliers = detect_outliers(x)
    samp.dat.outlier[i[outliers]] .= true
end

"""
    moving_median_indices(n::Int; b=2)

Generate index matrix for moving median calculation.

# Arguments
- `n`: Length of the time series
- `b`: Bandwidth (default: 2)

# Returns
- Matrix of indices where each row contains the indices for computing the moving median at that position
"""
function moving_median_indices(n::Int;b::Int=2)
    offset = [collect(-b:-1);collect(1:b)]
    i = collect(1:n)
    out = i .+ offset[1]
    for o in offset[2:end]
        out = hcat(out,i .+ o)
    end
    for j in 1:b
        out[j,1:b+1-j] .= out[j,end] .+ 1 .- out[j,1:b+1-j]
    end
    for j in n-b+1:n
        out[j,b+1+n-j:2*b] .= out[j,1] .- collect(1:b+j-n)
    end
    return out
end
export moving_median_indices

function moving_median(vec::AbstractVector, i::Matrix{Int})
    n = length(vec)
    return [Statistics.median(vec[i[j, :]]) for j in 1:n]
end
function moving_median(mat::Matrix, i::Matrix{Int})
    out = copy(mat)
    for j in 1:size(mat,2)
        out[:,j] = moving_median(mat[:,j],i)
    end
    return out
end
