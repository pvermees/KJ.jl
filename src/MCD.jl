# Written by Gemini using the prompt
# "Implement the FAST-MCD algorithm in Julia."
# Minor subsequent edits by PV

"""
calculate_subset_stats(X, indices)

Calculates the mean, covariance matrix, and its 
determinant for a given subset of data.
"""
function calculate_subset_stats(X::AbstractMatrix, indices::Vector{Int})
    # Extract the subset
    X_H = X[indices, :]
    
    # Calculate Location (Mean) and Scatter (Covariance)
    mu = vec(Statistics.mean(X_H, dims=1))
    Sigma = Statistics.cov(X_H)
    
    # Calculate Determinant
    local d
    try
        d = LinearAlgebra.det(Sigma)
        if d <= 0
            d = Inf # Treat non-positive definite matrix as failed/singular
        end
    catch
        d = Inf # Handle singular matrix case
    end
    
    return mu, Sigma, d
end

"""
concentration_step(X, h, current_mu, current_Sigma)

The core C-Step. It finds the h observations closest to the current center (mu)
based on the Mahalanobis distance derived from the current scatter matrix (Sigma).
"""
function concentration_step(X::AbstractMatrix, h::Int, current_mu, current_Sigma)
    n, p = size(X)
    
    # If the current scatter matrix is singular or non-positive, we cannot proceed.
    if isinf(LinearAlgebra.det(current_Sigma)) ||
        LinearAlgebra.det(current_Sigma) < 1e-10
        return Int[], zeros(p), Matrix{eltype(X)}(I, p, p), Inf # Signal failure
    end

    # Calculate inverse of the scatter matrix
    Sigma_inv = inv(current_Sigma)
    
    # Calculate squared Mahalanobis distance for ALL n points
    distances_sq = zeros(n)
    for i in 1:n
        x_i = X[i, :]
        diff = x_i - current_mu
        # MD^2 = diff' * Sigma_inv * diff
        distances_sq[i] = dot(diff, Sigma_inv * diff)
    end
    
    # Select the h smallest distances (i.e., the h closest points)
    sorted_indices = sortperm(distances_sq)
    new_indices = sorted_indices[1:h]
    
    # Calculate the new location, scatter, and determinant from the new, tighter set
    new_mu, new_Sigma, new_det = calculate_subset_stats(X, new_indices)
    
    return new_indices, new_mu, new_Sigma, new_det
end

# --- 2. Main Algorithm ---

"""
fast_mcd(X::AbstractMatrix; h::Int=0, num_starts::Int=500, max_c_steps::Int=10)

Implements the FAST-MCD algorithm using randomized starts and concentration steps.

# Arguments
- `X`: The data matrix (observations as rows, features as columns).
- `h`: The size of the subset. Defaults to floor((n + p + 1) / 2) for max robustness.
- `num_starts`: Number of random initial subsets to check.
- `max_c_steps`: Maximum number of C-steps per start before declaring convergence.
"""
function fast_mcd(X::AbstractMatrix{T}; h::Int=0, num_starts::Int=500, max_c_steps::Int=10) where T <: Real
    n, p = size(X)

    # Set default 'h' for 50% breakdown point (maximum robustness)
    if h == 0
        h = floor(Int, (n + p + 1) / 2)
    end
    
    if h <= p
        error("Subset size h=$h must be greater than the number of features p=$p for non-singular covariance.")
    end

    # Initialize best results found across all starts
    min_det = Inf
    best_indices = Int[]
    best_location = zeros(T, p)
    best_scatter_raw = zeros(T, p, p)

    for i in 1:num_starts
        # 1. Initialization: Randomly select an initial subset of size h
        current_indices = Distributions.sample(1:n, h, replace=false)
        current_mu, current_Sigma, current_det =
            calculate_subset_stats(X, current_indices)
        
        # 2. Apply Concentration Steps (C-steps)
        c_step_count = 0
        while c_step_count < max_c_steps
            
            if isinf(current_det)
                break 
            end

            # Perform C-Step to find a tighter subset
            next_indices, next_mu, next_Sigma, next_det = concentration_step(
                X, h, current_mu, current_Sigma
            )

            # Check for convergence: if the determinant did not decrease, stop.
            if next_det < current_det
                # Improvement found: update and continue
                current_indices = next_indices
                current_mu = next_mu
                current_Sigma = next_Sigma
                current_det = next_det
                c_step_count += 1
            else
                # Convergence reached for this initial start
                break
            end
        end

        # 3. Final Selection: Update the overall best result found so far
        if current_det < min_det
            min_det = current_det
            best_indices = current_indices
            best_location = current_mu
            best_scatter_raw = current_Sigma
        end
    end
    
    if isempty(best_indices)
        error("FAST-MCD failed to find a valid non-singular subset " *
              "after $num_starts attempts.")
    end

    return (
        location = best_location,
        scatter_raw = best_scatter_raw,
        min_determinant = min_det,
        support_indices = best_indices
    )
end
export fast_mcd
