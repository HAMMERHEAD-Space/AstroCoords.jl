export Keplerian
"""
    Keplerian{T} <: AstroCoord

Keplerian Orbital Elements. 6D parameterziation of the orbit.

# Fields
a - semi-major axis
e - eccentricity
i - inclination
Ω - Right Ascension of Ascending Node
ω - Argument of Perigee
f - True Anomaly

# Exotic Properties
M - Mean Anomaly (computed from true anomaly and eccentricity)
E - Eccentric Anomaly (computed from true anomaly and eccentricity)

# Constructors
- `Keplerian(a, e, i, Ω, ω, f)`
- `Keplerian(X::AbstractArray)`
- `Keplerian(X::AstroCoord, μ::Number)`
"""
struct Keplerian{T} <: AstroCoord{6,T}
    a::T
    e::T
    i::T
    Ω::T
    ω::T
    f::T
    @inline Keplerian{T}(a, e, i, Ω, ω, f) where {T} = new{T}(a, e, i, Ω, ω, f)
    @inline Keplerian{T}(p::Keplerian{T}) where {T} = new{T}(p.a, p.e, p.i, p.Ω, p.ω, p.f)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
Keplerian(X::AbstractVector{T}) where {T} = Keplerian{T}(X[1], X[2], X[3], X[4], X[5], X[6])
function Keplerian(a::A, e::E, i::I, Ω::O, ω::W, f::F) where {A,E,I,O,W,F}
    return Keplerian{promote_type(A, E, I, O, W, F)}(a, e, i, Ω, ω, f)
end
# More specific than AbstractVector to avoid ambiguity
function Keplerian(g::StaticVector{N,T}) where {N,T}
    Keplerian{T}(g[1], g[2], g[3], g[4], g[5], g[6])
end
Keplerian{T}(g::StaticVector) where {T} = Keplerian{T}(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::Keplerian{T}) where {T<:Number} = SVector{6,T}(g.a, g.e, g.i, g.Ω, g.ω, g.f)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{K}; T::DataType=Float64) where {K<:Keplerian}
    return K{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Keplerian{T}, i::Int) where {T<:Number}
    if i == 1
        return p.a
    elseif i == 2
        return p.e
    elseif i == 3
        return p.i
    elseif i == 4
        return p.Ω
    elseif i == 5
        return p.ω
    elseif i == 6
        return p.f
    else
        throw(BoundsError(p, i))
    end
end

# ~~~~~~~~~~~~~~~ Property Interface for Exotic Coordinates ~~~~~~~~~~~~~~~ #
function Base.getproperty(kep::Keplerian, sym::Symbol)
    if sym === :M
        # Mean anomaly - computed from true anomaly and eccentricity
        return trueAnomaly2MeanAnomaly(kep.f, kep.e)
    elseif sym === :E
        # Eccentric anomaly - computed from true anomaly and eccentricity
        return trueAnomaly2EccentricAnomaly(kep.f, kep.e)
    else
        # Default behavior for regular fields
        return getfield(kep, sym)
    end
end

function Base.propertynames(::Keplerian)
    return (:a, :e, :i, :Ω, :ω, :f, :M, :E)
end
