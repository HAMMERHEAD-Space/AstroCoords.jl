# Core uncertainty propagation utilities
# Not imported by default - users activate this environment separately

using AstroCoords
using LinearAlgebra
using ForwardDiff
using Statistics
using Distributions
using StaticArrays

export propagate_linear, propagate_unscented, propagate_monte_carlo
export generate_sigma_points, reconstruct_statistics

"""
    propagate_linear(uc::UncertainCoord, transform_func, μ; atol=1e-10)

Propagate uncertainty through a coordinate transformation using first-order
linear approximation (Jacobian-based covariance propagation).

# Arguments
- `uc::UncertainCoord`: Input coordinate with uncertainty
- `transform_func`: Function that transforms coordinates, signature: f(coord::AstroCoord, μ) -> AstroCoord
- `μ::Number`: Gravitational parameter
- `atol::Number`: Absolute tolerance for numerical Jacobian computation

# Returns
- `UncertainCoord`: Transformed coordinate with propagated uncertainty

# Theory

For a nonlinear transformation y = f(x), the output covariance is approximated as:

Σ_y ≈ J * Σ_x * J^T

where J = ∂f/∂x is the Jacobian matrix. This is the first-order Taylor expansion
used in Extended Kalman Filters and standard orbit determination.

# Reference

Tapley, B. D., Schutz, B. E., & Born, G. H. (2004). *Statistical Orbit Determination*.
Academic Press. Chapter 4: "Estimation of Parameters".

# Example

```julia
using AstroCoords, LinearAlgebra

# Keplerian orbit with uncertainty
kep = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.5)
Σ_kep = Diagonal([100.0, 0.001, 1e-6, 1e-6, 1e-6, 1e-6])
unc_kep = UncertainCoord(kep, CovarianceUncertainty(Σ_kep))

# Propagate to Cartesian
μ = 398600.4418
transform = (k, μ) -> Cartesian(k, μ)
unc_cart = propagate_linear(unc_kep, transform, μ)
```
"""
function propagate_linear(uc::UncertainCoord, transform_func, μ; atol=1e-10)
    coord_in = uc.coord
    cov_in = uc.uncertainty.cov
    
    # Compute Jacobian using ForwardDiff
    x_in = params(coord_in)
    
    # Wrapper for ForwardDiff (needs to work with arrays)
    f_vec(x) = params(transform_func(typeof(coord_in)(x), μ))
    
    # Compute Jacobian
    J = ForwardDiff.jacobian(f_vec, x_in)
    
    # Propagate covariance: Σ_out = J * Σ_in * J^T
    cov_out = J * cov_in * J'
    
    # Ensure symmetry (numerical errors can break this)
    cov_out = 0.5 * (cov_out + cov_out')
    
    # Transform nominal coordinate
    coord_out = transform_func(coord_in, μ)
    
    return UncertainCoord(coord_out, CovarianceUncertainty(cov_out))
end

"""
    generate_sigma_points(mean::AbstractVector, cov::Matrix; α=1e-3, β=2.0, κ=0.0)

Generate sigma points for the Unscented Transform.

# Arguments
- `mean::AbstractVector`: Mean of the distribution (n-dimensional)
- `cov::Matrix`: Covariance matrix (n × n)
- `α::Number`: Spread parameter (typical: 1e-4 to 1)
- `β::Number`: Distribution parameter (β=2 optimal for Gaussian)
- `κ::Number`: Secondary scaling parameter (typically 0 or 3-n)

# Returns
- `sigma_points::Matrix`: 2n+1 sigma points (each column is one point)
- `weights_mean::Vector`: Weights for computing mean
- `weights_cov::Vector`: Weights for computing covariance

# Theory

The Unscented Transform approximates a probability distribution using 2n+1
carefully chosen sigma points:

χ₀ = μ
χᵢ = μ + (√((n+λ)Σ))ᵢ,  i = 1,...,n
χᵢ = μ - (√((n+λ)Σ))ᵢ₋ₙ,  i = n+1,...,2n

where λ = α²(n+κ) - n is a scaling parameter.

# Reference

Julier, S. J., & Uhlmann, J. K. (1997). "New extension of the Kalman filter to
nonlinear systems." *Signal Processing, Sensor Fusion, and Target Recognition VI*,
3068, 182-193.

# See also
- [`propagate_unscented`](@ref)
- [`reconstruct_statistics`](@ref)
"""
function generate_sigma_points(mean::AbstractVector, cov::Matrix; α=1e-3, β=2.0, κ=0.0)
    n = length(mean)
    λ = α^2 * (n + κ) - n
    
    # Compute matrix square root: √((n+λ)Σ)
    # Use Cholesky decomposition for numerical stability
    sqrt_factor = √(n + λ)
    L = cholesky(Hermitian(cov)).L
    sqrt_cov = sqrt_factor * L
    
    # Initialize sigma points
    sigma_points = zeros(n, 2n + 1)
    
    # First point is the mean
    sigma_points[:, 1] = mean
    
    # Next n points: mean + columns of sqrt_cov
    for i in 1:n
        sigma_points[:, i+1] = mean + sqrt_cov[:, i]
    end
    
    # Last n points: mean - columns of sqrt_cov
    for i in 1:n
        sigma_points[:, n+i+1] = mean - sqrt_cov[:, i]
    end
    
    # Compute weights
    weights_mean = zeros(2n + 1)
    weights_cov = zeros(2n + 1)
    
    weights_mean[1] = λ / (n + λ)
    weights_cov[1] = λ / (n + λ) + (1 - α^2 + β)
    
    for i in 2:(2n+1)
        weights_mean[i] = 1 / (2 * (n + λ))
        weights_cov[i] = 1 / (2 * (n + λ))
    end
    
    return sigma_points, weights_mean, weights_cov
end

"""
    reconstruct_statistics(transformed_points::Matrix, weights_mean::Vector, weights_cov::Vector)

Reconstruct mean and covariance from transformed sigma points.

# Arguments
- `transformed_points::Matrix`: Transformed sigma points (n × (2n+1))
- `weights_mean::Vector`: Weights for mean computation
- `weights_cov::Vector`: Weights for covariance computation

# Returns
- `mean::Vector`: Reconstructed mean
- `cov::Matrix`: Reconstructed covariance

# See also
- [`generate_sigma_points`](@ref)
- [`propagate_unscented`](@ref)
"""
function reconstruct_statistics(transformed_points::Matrix, weights_mean::Vector, weights_cov::Vector)
    n, num_points = size(transformed_points)
    
    # Compute weighted mean
    mean = zeros(n)
    for i in 1:num_points
        mean += weights_mean[i] * transformed_points[:, i]
    end
    
    # Compute weighted covariance
    cov = zeros(n, n)
    for i in 1:num_points
        diff = transformed_points[:, i] - mean
        cov += weights_cov[i] * (diff * diff')
    end
    
    # Ensure symmetry
    cov = 0.5 * (cov + cov')
    
    return mean, cov
end

"""
    propagate_unscented(uc::UncertainCoord, transform_func, μ; α=1e-3, β=2.0, κ=0.0)

Propagate uncertainty through a coordinate transformation using the Unscented Transform.

# Arguments
- `uc::UncertainCoord`: Input coordinate with uncertainty
- `transform_func`: Transformation function
- `μ::Number`: Gravitational parameter
- `α, β, κ`: Unscented Transform tuning parameters

# Returns
- `UncertainCoord`: Transformed coordinate with propagated uncertainty

# Theory

The Unscented Transform propagates 2n+1 sigma points through the nonlinear
transformation, then reconstructs the output statistics. This captures nonlinearity
more accurately than linear approximation, with computational cost O(n³) vs O(n²)
for gradient-based methods.

# Reference

Julier, S. J., & Uhlmann, J. K. (1997). *Signal Processing, Sensor Fusion, and
Target Recognition VI*, 3068, 182-193.

# Example

```julia
unc_kep = UncertainCoord(kep, CovarianceUncertainty(Σ_kep))
transform = (k, μ) -> Cartesian(k, μ)
unc_cart = propagate_unscented(unc_kep, transform, μ)
```
"""
function propagate_unscented(uc::UncertainCoord, transform_func, μ; α=1e-3, β=2.0, κ=0.0)
    coord_in = uc.coord
    cov_in = uc.uncertainty.cov
    
    # Generate sigma points
    x_in = Vector(params(coord_in))
    sigma_points, weights_mean, weights_cov = generate_sigma_points(x_in, cov_in; α=α, β=β, κ=κ)
    
    # Propagate sigma points through transformation
    n = length(x_in)
    num_points = 2n + 1
    
    # Transform first point to get output dimension
    test_out = transform_func(typeof(coord_in)(sigma_points[:, 1]), μ)
    m = length(params(test_out))
    
    transformed_points = zeros(m, num_points)
    
    for i in 1:num_points
        coord_i = typeof(coord_in)(sigma_points[:, i])
        transformed_i = transform_func(coord_i, μ)
        transformed_points[:, i] = params(transformed_i)
    end
    
    # Reconstruct statistics
    mean_out, cov_out = reconstruct_statistics(transformed_points, weights_mean, weights_cov)
    
    # Create output coordinate from mean
    CoordType = typeof(test_out)
    coord_out = CoordType(mean_out)
    
    return UncertainCoord(coord_out, CovarianceUncertainty(cov_out))
end

"""
    propagate_monte_carlo(uc::UncertainCoord, transform_func, μ; n_samples=10000, rng=Random.GLOBAL_RNG)

Propagate uncertainty using Monte Carlo sampling.

# Arguments
- `uc::UncertainCoord`: Input coordinate with uncertainty
- `transform_func`: Transformation function
- `μ::Number`: Gravitational parameter
- `n_samples::Int`: Number of Monte Carlo samples
- `rng`: Random number generator

# Returns
- `UncertainCoord`: Transformed coordinate with propagated uncertainty
- `samples::Matrix`: Matrix of output samples (m × n_samples)

# Theory

Monte Carlo propagation:
1. Sample from input distribution: xᵢ ~ N(μ, Σ)
2. Transform each sample: yᵢ = f(xᵢ)
3. Estimate output statistics from {yᵢ}

Convergence rate: O(1/√N), so error decreases slowly with sample count.

# Example

```julia
unc_kep = UncertainCoord(kep, CovarianceUncertainty(Σ_kep))
transform = (k, μ) -> Cartesian(k, μ)
unc_cart, samples = propagate_monte_carlo(unc_kep, transform, μ; n_samples=10000)
```
"""
function propagate_monte_carlo(uc::UncertainCoord, transform_func, μ; n_samples=10000, rng=Random.GLOBAL_RNG)
    coord_in = uc.coord
    cov_in = uc.uncertainty.cov
    
    # Create multivariate normal distribution
    x_in = Vector(params(coord_in))
    dist = MvNormal(x_in, Hermitian(cov_in))
    
    # Sample from input distribution
    samples_in = rand(rng, dist, n_samples)
    
    # Transform first sample to get output dimension
    test_out = transform_func(typeof(coord_in)(samples_in[:, 1]), μ)
    m = length(params(test_out))
    
    # Transform all samples
    samples_out = zeros(m, n_samples)
    
    for i in 1:n_samples
        coord_i = typeof(coord_in)(samples_in[:, i])
        transformed_i = transform_func(coord_i, μ)
        samples_out[:, i] = params(transformed_i)
    end
    
    # Estimate statistics from samples
    mean_out = vec(mean(samples_out, dims=2))
    cov_out = cov(samples_out, dims=2)
    
    # Create output coordinate
    CoordType = typeof(test_out)
    coord_out = CoordType(mean_out)
    
    return UncertainCoord(coord_out, CovarianceUncertainty(cov_out)), samples_out
end