export EDromo
"""
    EDromo{T} <: AstroCoord

EDromo state vector.
ζ₁ - In-plane element 1
ζ₂ - In-plane element 2
ζ₃ - Corresponds to orbital energy
ζ₄ - Quaternion element 1
ζ₅ - Quaternion element 2
ζ₆ - Quaternion element 3
ζ₇ - Quaternion element 4
ζ₈ - Time element

Constructors
EDromo(ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈)
EDromo(X::AbstractArray)

"""
struct EDromo{T} <: AstroCoord{8,T}
    ζ₁::T
    ζ₂::T
    ζ₃::T
    ζ₄::T
    ζ₅::T
    ζ₆::T
    ζ₇::T
    ζ₈::T
    @inline function EDromo{T}(ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈) where {T}
        return new{T}(ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈)
    end
    @inline EDromo{T}(p::EDromo{T}) where {T} = new{T}(
        p.ζ₁, p.ζ₂, p.ζ₃, p.ζ₄, p.ζ₅, p.ζ₆, p.ζ₇, p.ζ₈
    )
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
function EDromo(X::AbstractVector{T}) where {T}
    EDromo{T}(X[1], X[2], X[3], X[4], X[5], X[6], X[7], X[8])
end
function EDromo(
    ζ₁::T1, ζ₂::T2, ζ₃::T3, ζ₄::T4, ζ₅::T5, ζ₆::T6, ζ₇::T7, ζ₈::T8
) where {T1,T2,T3,T4,T5,T6,T7,T8}
    T = promote_type(T1, T2, T3, T4, T5, T6, T7, T8)
    return EDromo{T}(T(ζ₁), T(ζ₂), T(ζ₃), T(ζ₄), T(ζ₅), T(ζ₆), T(ζ₇), T(ζ₈))
end
function (::Type{E})(g::StaticVector) where {E<:EDromo}
    E(g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8])
end

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
function params(g::EDromo{T}) where {T<:Number}
    SVector{8,T}(g.ζ₁, g.ζ₂, g.ζ₃, g.ζ₄, g.ζ₅, g.ζ₆, g.ζ₇, g.ζ₈)
end

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{E}; T::DataType=Float64) where {E<:EDromo}
    return E{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::EDromo{T}, i::Int) where {T<:Number}
    if i == 1
        return p.ζ₁
    elseif i == 2
        return p.ζ₂
    elseif i == 3
        return p.ζ₃
    elseif i == 4
        return p.ζ₄
    elseif i == 5
        return p.ζ₅
    elseif i == 6
        return p.ζ₆
    elseif i == 7
        return p.ζ₇
    elseif i == 8
        return p.ζ₈
    else
        throw(BoundsError(p, i))
    end
end
