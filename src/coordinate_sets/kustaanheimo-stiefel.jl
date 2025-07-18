export KustaanheimoStiefel
"""
    KustaanheimoStiefel{T} <: AstroCoord

Kustaanheimo-Stiefel (KS) state vector. The KS transformation regularizes the 
two-body problem by transforming the 3D Cartesian coordinates into a 4D space, 
which eliminates the singularity at `r=0`.

This state vector contains 10 elements.
u₁ - First component of the KS position vector
u₂ - Second component of the KS position vector
u₃ - Third component of the KS position vector
u₄ - Fourth component of the KS position vector
u₁_prime - First component of the KS velocity vector
u₂_prime - Second component of the KS velocity vector
u₃_prime - Third component of the KS velocity vector
u₄_prime - Fourth component of the KS velocity vector
h - The negative of the total orbital energy (`-E`)
τ - The time element, which can be either physical time `t` or a linear time element.

Constructors
KustaanheimoStiefel(u₁, u₂, u₃, u₄, u₁_prime, u₂_prime, u₃_prime, u₄_prime, h, τ)
KustaanheimoStiefel(X::AbstractArray)

"""
struct KustaanheimoStiefel{T} <: AstroCoord{10,T}
    u₁::T
    u₂::T
    u₃::T
    u₄::T
    u₁_prime::T
    u₂_prime::T
    u₃_prime::T
    u₄_prime::T
    h::T
    τ::T
    @inline function KustaanheimoStiefel{T}(
        u₁, u₂, u₃, u₄, u₁_prime, u₂_prime, u₃_prime, u₄_prime, h, τ
    ) where {T}
        return new{T}(u₁, u₂, u₃, u₄, u₁_prime, u₂_prime, u₃_prime, u₄_prime, h, τ)
    end
    @inline KustaanheimoStiefel{T}(p::KustaanheimoStiefel{T}) where {T} = new{T}(
        p.u₁, p.u₂, p.u₃, p.u₄, p.u₁_prime, p.u₂_prime, p.u₃_prime, p.u₄_prime, p.h, p.τ
    )
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
function KustaanheimoStiefel(X::AbstractVector{T}) where {T}
    KustaanheimoStiefel{T}(X[1], X[2], X[3], X[4], X[5], X[6], X[7], X[8], X[9], X[10])
end

function KustaanheimoStiefel(
    u₁::T1,
    u₂::T2,
    u₃::T3,
    u₄::T4,
    u₁_prime::T5,
    u₂_prime::T6,
    u₃_prime::T7,
    u₄_prime::T8,
    h::T9,
    τ::T10,
) where {T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
    T = promote_type(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10)
    return KustaanheimoStiefel{T}(
        u₁, u₂, u₃, u₄, u₁_prime, u₂_prime, u₃_prime, u₄_prime, h, τ
    )
end
function (::Type{KS})(g::StaticVector) where {KS<:KustaanheimoStiefel}
    KS(g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9], g[10])
end

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
function params(g::KustaanheimoStiefel{T}) where {T<:Number}
    SVector{10,T}(
        g.u₁, g.u₂, g.u₃, g.u₄, g.u₁_prime, g.u₂_prime, g.u₃_prime, g.u₄_prime, g.h, g.τ
    )
end

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{KS}; T::DataType=Float64) where {KS<:KustaanheimoStiefel}
    return KS{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::KustaanheimoStiefel{T}, i::Int) where {T<:Number}
    if i == 1
        return p.u₁
    elseif i == 2
        return p.u₂
    elseif i == 3
        return p.u₃
    elseif i == 4
        return p.u₄
    elseif i == 5
        return p.u₁_prime
    elseif i == 6
        return p.u₂_prime
    elseif i == 7
        return p.u₃_prime
    elseif i == 8
        return p.u₄_prime
    elseif i == 9
        return p.h
    elseif i == 10
        return p.τ
    else
        throw(BoundsError(p, i))
    end
end
