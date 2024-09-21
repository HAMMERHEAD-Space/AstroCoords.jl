export Delauney
"""
    Delauney{T} <: AstroCoord

Delauney Orbital Elements. 6D parameterziation of the orbit.
L - Canonical Keplerian Energy
G - Canonical Total Angular Momentum
H - Canonical Normal Angular Momentum (Relative to Equator)
M - Mean Anomaly
ω - Arugment of Periapsis
Ω - Right Ascention of the Ascending Node

Constructors
Delauney(L, G, H, M, ω, Ω)
Delauney(X::AbstractVector{<:Number})
Delauney(X::AstroCoord, μ::Number)

"""
struct Delauney{T} <: AstroCoord{6,T}
    L::T
    G::T
    H::T
    M::T
    ω::T
    Ω::T
    @inline Delauney{T}(L, G, H, M, ω, Ω) where {T} = new{T}(L, G, H, M, ω, Ω)
    @inline Delauney{T}(p::Delauney) where {T} = new{T}(p.L, p.G, p.H, p.M, p.ω, p.Ω)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
Delauney(X::AbstractVector{T}) where {T} = Cartesian{T}(X...)
function Delauney(L::LT, G::GT, H::HT, M::MT, ω::OT, Ω::OmT) where {LT,GT,HT,MT,OT,OmT}
    return Delauney{promote_type(LT, GT, HT, MT, OT, Omt)}(L, G, H, M, ω, Ω)
end
(::Type{D})(g::StaticVector) where {D<:Delauney} = D(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::Delauney) = SVector{6}(g.L, g.G, g.H, g.M, g.ω, g.Ω)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.one(::Type{D}) where {D<:Delauney} = D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Delauney, i::Int)
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
        throw(BoundsError(r, i))
    end
end
