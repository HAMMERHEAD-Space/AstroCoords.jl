export StiefelScheifele
"""
    StiefelScheifele{T} <: AstroCoord

Stiefel-Scheifele orbital elements. This is a 10-dimensional parameterization
of a regularized orbit.

The elements are:
α1 - The first component of the Stiefel-Scheifele position vector
α2 - The second component of the Stiefel-Scheifele position vector
α3 - The third component of the Stiefel-Scheifele position vector
α4 - The fourth component of the Stiefel-Scheifele position vector
β1 - The first component of the Stiefel-Scheifele velocity vector
β2 - The second component of the Stiefel-Scheifele velocity vector
β3 - The third component of the Stiefel-Scheifele velocity vector
β4 - The fourth component of the Stiefel-Scheifele velocity vector
ω - Related to the total energy of the orbit
t - Time or a time element, depending on `flag_time`

Constructors:
    StiefelScheifele(α1, α2, α3, α4, β1, β2, β3, β4, ω, t)
    StiefelScheifele(v::AbstractVector)
"""
struct StiefelScheifele{T} <: AstroCoord{10,T}
    α1::T
    α2::T
    α3::T
    α4::T
    β1::T
    β2::T
    β3::T
    β4::T
    ω::T
    t::T

    @inline function StiefelScheifele{T}(α1, α2, α3, α4, β1, β2, β3, β4, ω, t) where {T}
        return new{T}(α1, α2, α3, α4, β1, β2, β3, β4, ω, t)
    end
    @inline function StiefelScheifele{T}(v::StiefelScheifele{T}) where {T}
        return new{T}(v.α1, v.α2, v.α3, v.α4, v.β1, v.β2, v.β3, v.β4, v.ω, v.t)
    end
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
function StiefelScheifele(X::AbstractVector{T}) where {T}
    StiefelScheifele{T}(X[1], X[2], X[3], X[4], X[5], X[6], X[7], X[8], X[9], X[10])
end
function StiefelScheifele(
    α1::A1, α2::A2, α3::A3, α4::A4, β1::B1, β2::B2, β3::B3, β4::B4, ω::W, t::TT
) where {A1,A2,A3,A4,B1,B2,B3,B4,W,TT}
    T = promote_type(A1, A2, A3, A4, B1, B2, B3, B4, W, TT)
    return StiefelScheifele{T}(α1, α2, α3, α4, β1, β2, β3, β4, ω, t)
end
function (::Type{SS})(g::StaticVector) where {SS<:StiefelScheifele}
    SS(g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9], g[10])
end

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
function params(g::StiefelScheifele{T}) where {T<:Number}
    return SVector{10,T}(g.α1, g.α2, g.α3, g.α4, g.β1, g.β2, g.β3, g.β4, g.ω, g.t)
end

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{SS}; T::DataType=Float64) where {SS<:StiefelScheifele}
    return SS{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::StiefelScheifele, i::Int)
    if i == 1
        return p.α1
    elseif i == 2
        return p.α2
    elseif i == 3
        return p.α3
    elseif i == 4
        return p.α4
    elseif i == 5
        return p.β1
    elseif i == 6
        return p.β2
    elseif i == 7
        return p.β3
    elseif i == 8
        return p.β4
    elseif i == 9
        return p.ω
    elseif i == 10
        return p.t
    else
        throw(BoundsError(p, i))
    end
end
