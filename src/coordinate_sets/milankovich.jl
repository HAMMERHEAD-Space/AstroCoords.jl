export Milankovich
"""
    Milankovich{T} <: AstroCoord

Milankovich Orbital Elements. 7D parameterziation of the orbit.
hx - X-component of Angular Momentum Vector
hy - Y-component of Angular Momentum Vector
hz - Z-component of Angular Momentum Vector
ex - X-component of Eccentricity Vector
ey - Y-component of Eccentricity Vector
ez - Z-component of Eccentricity Vector
L - True Longitude

Constructors
Milankovich(hx, hy, hz, ex, ey, ez, L)
Milankovich(X::AbstractArray)
Milankovich(X::AstroCoord, μ::Number)

"""
struct Milankovich{T} <: AstroCoord{7,T}
    hx::T
    hy::T
    hz::T
    ex::T
    ey::T
    ez::T
    L::T
    @inline Milankovich{T}(hx, hy, hz, ex, ey, ez, L) where {T} = new{T}(
        hx, hy, hz, ex, ey, ez, L
    )
    @inline Milankovich{T}(X::Milankovich{T}) where {T} = new{T}(
        X.hx, X.hy, X.hz, X.ex, X.ey, X.ez, X.L
    )
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
function Milankovich(X::AbstractVector{T}) where {T}
    return Milankovich{T}(X[1], X[2], X[3], X[4], X[5], X[6], X[7])
end
function Milankovich(
    hx::HX, hy::HY, hz::HZ, ex::EX, ey::EY, ez::EZ, L::LT
) where {HX,HY,HZ,EX,EY,EZ,LT}
    return Milankovich{promote_type(HX, HY, HZ, EX, EY, EZ, LT)}(hx, hy, hz, ex, ey, ez, L)
end
# More specific than AbstractVector to avoid ambiguity
function Milankovich(g::StaticVector{N,T}) where {N,T}
    Milankovich{T}(g[1], g[2], g[3], g[4], g[5], g[6], g[7])
end
function Milankovich{T}(g::StaticVector) where {T}
    Milankovich{T}(g[1], g[2], g[3], g[4], g[5], g[6], g[7])
end

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
function params(M::Milankovich{T}) where {T<:Number}
    return SVector{7,T}(M.hx, M.hy, M.hz, M.ex, M.ey, M.ez, M.L)
end

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{M}; T::DataType=Float64) where {M<:Milankovich}
    return M{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::Milankovich{T}, i::Int) where {T<:Number}
    if i == 1
        return p.hx
    elseif i == 2
        return p.hy
    elseif i == 3
        return p.hz
    elseif i == 4
        return p.ex
    elseif i == 5
        return p.ey
    elseif i == 6
        return p.ez
    elseif i == 7
        return p.L
    else
        throw(BoundsError(p, i))
    end
end
