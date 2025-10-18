export Delaunay
"""
    Delaunay{T} <: AstroCoord

Delaunay Orbital Elements. 6D parameterziation of the orbit.
L - Canonical Keplerian Energy
G - Canonical Total Angular Momentum
H - Canonical Normal Angular Momentum (Relative to Equator)
M - Mean Anomaly
ω - Argument of Periapsis
Ω - Right Ascention of the Ascending Node

Constructors
Delaunay(L, G, H, M, ω, Ω)
Delaunay(X::AbstractVector{<:Number})
Delaunay(X::AstroCoord, μ::Number)

"""
struct Delaunay{T} <: AstroCoord{6,T}
    L::T
    G::T
    H::T
    M::T
    ω::T
    Ω::T
    @inline Delaunay{T}(L, G, H, M, ω, Ω) where {T} = new{T}(L, G, H, M, ω, Ω)
    @inline Delaunay{T}(p::Delaunay{T}) where {T} = new{T}(p.L, p.G, p.H, p.M, p.ω, p.Ω)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
Delaunay(X::AbstractVector{T}) where {T} = Delaunay{T}(X[1], X[2], X[3], X[4], X[5], X[6])
function Delaunay(L::LT, G::GT, H::HT, M::MT, ω::PT, Ω::OmT) where {LT,GT,HT,MT,PT,OmT}
    return Delaunay{promote_type(LT, GT, HT, MT, PT, OmT)}(L, G, H, M, ω, Ω)
end
# More specific than AbstractVector to avoid ambiguity
Delaunay(g::StaticVector{N,T}) where {N,T} = Delaunay{T}(g[1], g[2], g[3], g[4], g[5], g[6])
Delaunay{T}(g::StaticVector) where {T} = Delaunay{T}(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::Delaunay{T}) where {T<:Number} = SVector{6,T}(g.L, g.G, g.H, g.M, g.ω, g.Ω)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{D}; T::DataType=Float64) where {D<:Delaunay}
    return D{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Delaunay{T}, i::Int) where {T<:Number}
    if i == 1
        return p.L
    elseif i == 2
        return p.G
    elseif i == 3
        return p.H
    elseif i == 4
        return p.M
    elseif i == 5
        return p.ω
    elseif i == 6
        return p.Ω
    else
        throw(BoundsError(p, i))
    end
end
