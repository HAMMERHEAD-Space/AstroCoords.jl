export Spherical
"""
    Spherical{T} <: AstroCoord

Spherical Orbital Elements. 6D parameterziation of the orbit
r - radius
θ - in-plane angle
ϕ - out-of-plane angle
ṙ - instantaneous rate of change of radius
θdot - instantaneous rate of change of θ
ϕdot - instantaneous rate of change of ϕ

Constructors
Spherical(r, θ, ϕ, ṙ, ω, Ω)
Spherical(X::AbstractArray)
Spherical(X::AstroCoord, μ::Number)

"""
struct Spherical{T} <: AstroCoord{6,T}
    r::T
    θ::T
    ϕ::T
    ṙ::T
    θdot::T
    ϕdot::T
    @inline Spherical{T}(r, θ, ϕ, ṙ, θdot, ϕdot) where {T} = new{T}(
        r, θ, ϕ, ṙ, θdot, ϕdot
    )
    @inline Spherical{T}(X::Spherical{T}) where {T} = new{T}(
        X.r, X.θ, X.ϕ, X.ṙ, X.θdot, X.ϕdot
    )
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
Spherical(X::AbstractVector{T}) where {T} = Spherical{T}(X[1], X[2], X[3], X[4], X[5], X[6])
function Spherical(r::R, θ::T, ϕ::P, ṙ::RD, θdot::TD, ϕdot::PD) where {R,T,P,RD,TD,PD}
    return Spherical{promote_type(R, T, P, RD, TD, PD)}(r, θ, ϕ, ṙ, θdot, ϕdot)
end
function (::Type{S})(g::StaticVector) where {S<:Spherical}
    return S(g[1], g[2], g[3], g[4], g[5], g[6])
end

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
function params(g::Spherical{T}) where {T<:Number}
    return SVector{6,T}(g.r, g.θ, g.ϕ, g.ṙ, g.θdot, g.ϕdot)
end

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{S}; T::DataType=Float64) where {S<:Spherical}
    return S{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Spherical{T}, i::Int) where {T<:Number}
    if i == 1
        return p.r
    elseif i == 2
        return p.θ
    elseif i == 3
        return p.ϕ
    elseif i == 4
        return p.ṙ
    elseif i == 5
        return p.θdot
    elseif i == 6
        return p.ϕdot
    else
        throw(BoundsError(p, i))
    end
end

export Cylindrical
"""
    Cylindrical{T} <: AstroCoord

Cylindrical Orbital Elements. 6D parameterziation of the orbit
ρ - in-plane radius
θ - in-plane angle
z - out-of-plane distance
ρdot - instantaneous rate of change of in-plsne radius
θdot - instantaneous rate of change of θ
ż - instantaneous rate of change of z

Constructors
Cylindrical(r, θ, z, ṙ, θdot, ż)
Cylindrical(X::AbstractArray)
Cylindrical(X::AstroCoord, μ::Number)

"""
struct Cylindrical{T} <: AstroCoord{6,T}
    ρ::T
    θ::T
    z::T
    ρdot::T
    θdot::T
    ż::T
    @inline Cylindrical{T}(ρ, θ, z, ρdot, θdot, ż) where {T} = new{T}(
        ρ, θ, z, ρdot, θdot, ż
    )
    @inline Cylindrical{T}(X::Cylindrical{T}) where {T} = new{T}(
        X.ρ, X.θ, X.z, X.ρdot, X.θdot, X.ż
    )
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
function Cylindrical(X::AbstractVector{T}) where {T}
    return Cylindrical{T}(X[1], X[2], X[3], X[4], X[5], X[6])
end
function Cylindrical(ρ::R, θ::T, z::Z, ρdot::RD, θdot::TD, ż::ZD) where {R,T,Z,RD,TD,ZD}
    return Cylindrical{promote_type(R, T, Z, RD, TD, ZD)}(ρ, θ, z, ρdot, θdot, ż)
end
function (::Type{C})(g::StaticVector) where {C<:Cylindrical}
    return C(g[1], g[2], g[3], g[4], g[5], g[6])
end

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
function params(g::Cylindrical{T}) where {T<:Number}
    return SVector{6,T}(g.ρ, g.θ, g.z, g.ρdot, g.θdot, g.ż)
end

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{C}; T::DataType=Float64) where {C<:Cylindrical}
    return C{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Cylindrical{T}, i::Int) where {T<:Number}
    if i == 1
        return p.ρ
    elseif i == 2
        return p.θ
    elseif i == 3
        return p.z
    elseif i == 4
        return p.ρdot
    elseif i == 5
        return p.θdot
    elseif i == 6
        return p.ż
    else
        throw(BoundsError(p, i))
    end
end
