export Cartesian
"""
    Cartesian{T} <: AstroCoord

Cartesian Orbital Elements. 6D parameterziation of the orbit.
x - X-position
y - Y-position
z - Z-position
ẋ - X-velocity
ẏ - Y-velocity
ż - Z-velocity

Constructors
Cartesian(x, y, z, ẋ, ẏ, ż)
Cartesian(X::AbstractArray)
Cartesian(X::AstroCoord, μ::Number)

"""
struct Cartesian{T} <: AstroCoord{6,T}
    x::T
    y::T
    z::T
    ẋ::T
    ẏ::T
    ż::T
    @inline Cartesian{T}(x, y, z, ẋ, ẏ, ż) where {T} = new{T}(x, y, z, ẋ, ẏ, ż)
    @inline Cartesian{T}(p::Cartesian{T}) where {T} = new{T}(
        p.x, p.y, p.z, p.ẋ, p.ẏ, p.ż
    )
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
Cartesian(X::AbstractVector{T}) where {T} = Cartesian{T}(X[1], X[2], X[3], X[4], X[5], X[6])
function Cartesian(x::X, y::Y, z::Z, ẋ::XV, ẏ::YV, ż::ZV) where {X,Y,Z,XV,YV,ZV}
    return Cartesian{promote_type(X, Y, Z, XV, YV, ZV)}(x, y, z, ẋ, ẏ, ż)
end
# More specific than AbstractVector to avoid ambiguity
function Cartesian(g::StaticVector{N,T}) where {N,T}
    Cartesian{T}(g[1], g[2], g[3], g[4], g[5], g[6])
end
Cartesian{T}(g::StaticVector) where {T} = Cartesian{T}(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::Cartesian{T}) where {T<:Number} = SVector{6,T}(g.x, g.y, g.z, g.ẋ, g.ẏ, g.ż)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{C}; T::DataType=Float64) where {C<:Cartesian}
    return C{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Cartesian{T}, i::Int) where {T<:Number}
    if i == 1
        return p.x
    elseif i == 2
        return p.y
    elseif i == 3
        return p.z
    elseif i == 4
        return p.ẋ
    elseif i == 5
        return p.ẏ
    elseif i == 6
        return p.ż
    else
        throw(BoundsError(p, i))
    end
end

# ~~~~~~~~~~~~~~~ Property Interface for Exotic Coordinates ~~~~~~~~~~~~~~~ #
function Base.getproperty(cart::Cartesian, sym::Symbol)
    if sym === :r
        return SVector{3}(cart.x, cart.y, cart.z)
    elseif sym === :v
        return SVector{3}(cart.ẋ, cart.ẏ, cart.ż)
    else
        # Default behavior for regular fields
        return getfield(cart, sym)
    end
end

function Base.propertynames(::Cartesian)
    return (:x, :y, :z, :ẋ, :ẏ, :ż, :r, :v)
end