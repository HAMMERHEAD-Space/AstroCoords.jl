export Poincare
"""
    Poincare{T} <: AstroCoord

Poincaré Canonical Orbital Elements. 6D parametrization of the orbit using 
Cartesian-style canonical coordinate-momentum pairs derived from Delaunay variables.

These elements are non-singular for circular orbits (e → 0) and equatorial orbits (i → 0),
making them particularly useful in perturbation theory and Hamiltonian celestial mechanics.

Supports both elliptic (a > 0) and hyperbolic (a < 0) orbits. For hyperbolic orbits,
the eccentricity action is defined as P = Λ + G (instead of Λ - G for elliptic) to keep 
all quantities real-valued. The orbit type is recovered unambiguously during inversion:
P ≤ Λ → elliptic, P > Λ → hyperbolic.

# Fields
Λ - Canonical action for semi-major axis: Λ = √(μ|a|)
λ - Mean longitude: λ = M + ω + Ω
ξ - Eccentricity cosine component: ξ = √(2P) cos(ω̃)
η - Eccentricity sine component: η = -√(2P) sin(ω̃), where ω̃ = ω + Ω
p - Inclination cosine component: p = √(2Q) cos(Ω), where Q = G(1 - cos(i))
q - Inclination sine component: q = -√(2Q) sin(Ω)

where P = Λ - G (elliptic) or P = Λ + G (hyperbolic), and G is the Delaunay angular momentum.

# Constructors
- `Poincare(Λ, λ, ξ, η, p, q)`
- `Poincare(X::AbstractVector{<:Number})`
- `Poincare(X::AstroCoord, μ::Number)`

# References
- Murray, C.D. and Dermott, S.F. "Solar System Dynamics." Cambridge University Press (1999).
- Laskar, J. and Robutel, P. "Stability of the Planetary Three-Body Problem." Celestial Mechanics and Dynamical Astronomy 62 (1995): 193-217.
"""
struct Poincare{T} <: AstroCoord{6,T}
    Λ::T
    λ::T
    ξ::T
    η::T
    p::T
    q::T
    @inline Poincare{T}(Λ, λ, ξ, η, p, q) where {T} = new{T}(Λ, λ, ξ, η, p, q)
    @inline Poincare{T}(X::Poincare{T}) where {T} = new{T}(X.Λ, X.λ, X.ξ, X.η, X.p, X.q)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
Poincare(X::AbstractVector{T}) where {T} = Poincare{T}(X[1], X[2], X[3], X[4], X[5], X[6])
function Poincare(Λ::LT, λ::LT2, ξ::XT, η::ET, p::PT, q::QT) where {LT,LT2,XT,ET,PT,QT}
    return Poincare{promote_type(LT, LT2, XT, ET, PT, QT)}(Λ, λ, ξ, η, p, q)
end
# More specific than AbstractVector to avoid ambiguity
Poincare(g::StaticVector{N,T}) where {N,T} = Poincare{T}(g[1], g[2], g[3], g[4], g[5], g[6])
Poincare{T}(g::StaticVector) where {T} = Poincare{T}(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::Poincare{T}) where {T<:Number} = SVector{6,T}(g.Λ, g.λ, g.ξ, g.η, g.p, g.q)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{P}; T::DataType=Float64) where {P<:Poincare}
    return P{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(pc::Poincare{T}, i::Int) where {T<:Number}
    if i == 1
        return pc.Λ
    elseif i == 2
        return pc.λ
    elseif i == 3
        return pc.ξ
    elseif i == 4
        return pc.η
    elseif i == 5
        return pc.p
    elseif i == 6
        return pc.q
    else
        throw(BoundsError(pc, i))
    end
end
