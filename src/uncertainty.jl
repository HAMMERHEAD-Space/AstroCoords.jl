export AbstractUncertainty, CovarianceUncertainty, UncertainCoord

"""
    abstract type AbstractUncertainty

Abstract base type for uncertainty representations. Subtypes should implement
methods for propagating uncertainty through coordinate transformations.

# Extended help

This type is designed to be extensible. Users can define custom uncertainty
representations by subtyping `AbstractUncertainty` and implementing appropriate
propagation methods.

## Examples of subtypes
- `CovarianceUncertainty`: Gaussian uncertainty represented by covariance matrix
- Custom distribution types from Distributions.jl
- Mixture models, particle representations, etc.

# See also
- [`CovarianceUncertainty`](@ref)
- [`UncertainCoord`](@ref)
"""
abstract type AbstractUncertainty end

"""
    CovarianceUncertainty{T} <: AbstractUncertainty

Represents Gaussian uncertainty via a covariance matrix.

# Fields
- `cov::Matrix{T}`: Covariance matrix (must be symmetric positive semi-definite)

# Examples
```julia
using LinearAlgebra

# 3x3 covariance for position uncertainty (km²)
Σ_pos = Diagonal([1.0, 1.0, 1.0])  # 1 km std in each direction
unc = CovarianceUncertainty(Σ_pos)
```

# See also
- [`AbstractUncertainty`](@ref)
- [`UncertainCoord`](@ref)
"""
struct CovarianceUncertainty{T<:Real} <: AbstractUncertainty
    cov::Matrix{T}

    function CovarianceUncertainty{T}(cov::Matrix{T}) where {T<:Real}
        if !issymmetric(cov)
            throw(ArgumentError("Covariance matrix must be symmetric"))
        end
        new{T}(cov)
    end
end

function CovarianceUncertainty(cov::Matrix{T}) where {T<:Real}
    return CovarianceUncertainty{T}(cov)
end

"""
    UncertainCoord{C<:AstroCoord, U<:AbstractUncertainty}

Wraps an astrodynamics coordinate with an uncertainty representation.

# Fields
- `coord::C`: The nominal coordinate value
- `uncertainty::U`: The uncertainty representation

# Examples
```julia
using LinearAlgebra

# Create a Keplerian orbit with uncertainty in semi-major axis and eccentricity
kep = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.5)  # km, rad
Σ = Diagonal([100.0, 0.001, 0.0, 0.0, 0.0, 0.0])  # 10 km std in a, 0.001 in e
unc_kep = UncertainCoord(kep, CovarianceUncertainty(Σ))

# Transform to Cartesian (uncertainty propagates)
μ = 398600.4418  # km³/s²
unc_cart = UncertainCoord(Cartesian(unc_kep.coord, μ), propagate_covariance(unc_kep, μ))
```

# Notes

This is a minimal implementation suitable for the core library. For advanced
uncertainty propagation methods (Unscented Transform, Monte Carlo, etc.),
see the examples in `examples/probabilistic_programming/`.

# See also
- [`AbstractUncertainty`](@ref)
- [`CovarianceUncertainty`](@ref)
"""
struct UncertainCoord{C<:AstroCoord,U<:AbstractUncertainty}
    coord::C
    uncertainty::U
end

# Accessor for coordinate
function Base.getproperty(uc::UncertainCoord, sym::Symbol)
    if sym === :coord || sym === :uncertainty
        getfield(uc, sym)
    else
        getproperty(getfield(uc, :coord), sym)
    end
end
