export GEqOE

"""
    GEqOE{T} <: AstroCoord

Generalized Equinoctial Orbital Elements (GEqOE). 6D parametrization of the orbit
that generalizes the equinoctial elements when perturbing forces are derived from
a disturbing potential.

Reference: Baù, G., Hernando-Ayuso, J., & Bombardelli, C. (2021).
"A generalization of the equinoctial orbital elements."
Celestial Mechanics and Dynamical Astronomy, 133(9), 1-32.

# Fields
- ν::T - Generalized mean motion: ν = (1/μ)(-2E)^(3/2)
- p₁::T - Generalized eccentricity component: p₁ = g sin(Ψ)
- p₂::T - Generalized eccentricity component: p₂ = g cos(Ψ)
- L::T - Generalized mean longitude
- q₁::T - Inclination component: q₁ = tan(i/2) sin(Ω)
- q₂::T - Inclination component: q₂ = tan(i/2) cos(Ω)

Where:
- E is the total energy (including disturbing potential U)
- g is the generalized Laplace vector magnitude
- Ψ is the generalized longitude of pericenter
- i is the orbital inclination
- Ω is the longitude of the ascending node

# Constructors
- `GEqOE(ν, p₁, p₂, L, q₁, q₂)`
- `GEqOE(X::AbstractArray)`
- `GEqOE(X::Cartesian, μ::Number, config::RegularizedCoordinateConfig)`

# Computed Properties
- `g`: Generalized eccentricity magnitude √(p₁² + p₂²)
- `Ψ`: Generalized longitude of pericenter atan(p₁, p₂)
- `L₀(geq, t)`: Generalized mean longitude at epoch L - νt

# Notes
- The perturbing potential `W` is passed via `RegularizedCoordinateConfig`. For Keplerian
  orbits use `W = 0`; for perturbed orbits precompute `W = V_total - V_keplerian` externally.
- Non-singular for circular and equatorial orbits
- Singular for retrograde equatorial orbits (i = π) and rectilinear motion
- Defined for negative total energy (E < 0)
"""
struct GEqOE{T} <: AstroCoord{6,T}
    ν::T
    p₁::T
    p₂::T
    L::T
    q₁::T
    q₂::T
    @inline GEqOE{T}(ν, p₁, p₂, L, q₁, q₂) where {T} = new{T}(ν, p₁, p₂, L, q₁, q₂)
    @inline GEqOE{T}(X::GEqOE{T}) where {T} = new{T}(X.ν, X.p₁, X.p₂, X.L, X.q₁, X.q₂)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
GEqOE(X::AbstractVector{T}) where {T} = GEqOE{T}(X[1], X[2], X[3], X[4], X[5], X[6])
function GEqOE(ν::NU, p₁::P1, p₂::P2, L::LT, q₁::Q1, q₂::Q2) where {NU,P1,P2,LT,Q1,Q2}
    return GEqOE{promote_type(NU, P1, P2, LT, Q1, Q2)}(ν, p₁, p₂, L, q₁, q₂)
end
# More specific than AbstractVector to avoid ambiguity
GEqOE(g::StaticVector{N,T}) where {N,T} = GEqOE{T}(g[1], g[2], g[3], g[4], g[5], g[6])
GEqOE{T}(g::StaticVector) where {T} = GEqOE{T}(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::GEqOE{T}) where {T<:Number} = SVector{6,T}(g.ν, g.p₁, g.p₂, g.L, g.q₁, g.q₂)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{G}; T::DataType=Float64) where {G<:GEqOE}
    return G{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::GEqOE{T}, i::Int) where {T<:Number}
    if i == 1
        return p.ν
    elseif i == 2
        return p.p₁
    elseif i == 3
        return p.p₂
    elseif i == 4
        return p.L
    elseif i == 5
        return p.q₁
    elseif i == 6
        return p.q₂
    else
        throw(BoundsError(p, i))
    end
end

# ~~~~~~~~~~~~~~~ Computed Properties ~~~~~~~~~~~~~~~ #
"""
Computed properties for GEqOE coordinate set.

# Additional Properties
- `a`: Generalized semi-major axis: a = (μ/ν²)^(1/3)
- `g`: Generalized eccentricity magnitude: g = √(p₁² + p₂²)
- `Ψ`: Generalized longitude of pericenter: Ψ = atan(p₁, p₂)
"""
function Base.getproperty(geq::GEqOE, sym::Symbol)
    if sym === :g
        # Generalized eccentricity magnitude: g = √(p₁² + p₂²)
        return sqrt(geq.p₁^2 + geq.p₂^2)
    elseif sym === :Ψ
        # Generalized longitude of pericenter: Ψ = atan(p₁, p₂)
        return atan(geq.p₁, geq.p₂)
    else
        return getfield(geq, sym)
    end
end

function Base.propertynames(::GEqOE)
    return (:ν, :p₁, :p₂, :L, :q₁, :q₂, :g, :Ψ)
end

# ~~~~~~~~~~~~~~~ Helper Functions ~~~~~~~~~~~~~~~ #
"""
    L₀(geq::GEqOE, t::Number; t₀::Number=0.0)

Generalized mean longitude at epoch: L₀ = L - νt
"""
L₀(geq::GEqOE, t::Number; t₀::Number=0.0) = geq.L - geq.ν * (t - t₀)

"""
    computed_a(geq::GEqOE, μ::Number)

Generalized semi-major axis: a = (μ/ν²)^(1/3)
"""
computed_a(geq::GEqOE, μ::Number) = (μ / geq.ν^2)^(1 / 3)

"""
    computed_ρ(geq::GEqOE, μ::Number)

Generalized semi-latus rectum: ρ = a(1 - g²)
"""
function computed_ρ(geq::GEqOE, μ::Number)
    a = computed_a(geq, μ)
    g² = geq.p₁^2 + geq.p₂^2
    return a * (1 - g²)
end
