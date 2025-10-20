export ModEq
"""
    ModEq{T} <: AstroCoord

Modified Equinoctial Orbital Elements. 6D parameterziation of the orbit.

# Fields
- `p` - semi-parameter
- `f` - eccentricity projection onto longitude of perigee
- `g` - eccentricity projection onto ⟂ longitude of perigee
- `h` - projection of half inclination onto RAAN
- `k` - projection of half inclination onto ⟂ RAAN
- `L` - true longitude

# Constructors
- `ModEq(p, f, g, h, k, L)`
- `ModEq(X::AbstractArray)`
- `ModEq(X::AstroCoord, μ::Number)`

# Examples
```julia
# Create ModEq coordinates
p = 6.93e6  # Semi-parameter
f = 0.1     # Eccentricity projection
modeq = ModEq(p, f, 0.0, 0.0, 0.0, π/3)
modeq.p  # 6.93e6 (semi-parameter)
modeq.f  # 0.1 (f component)
```

# See Also
- [`ModEqN`](@ref): Modified Equinoctial coordinates with mean motion instead of semi-parameter
"""
struct ModEq{T} <: AstroCoord{6,T}
    p::T
    f::T
    g::T
    h::T
    k::T
    L::T
    @inline ModEq{T}(p, f, g, h, k, L) where {T} = new{T}(p, f, g, h, k, L)
    @inline ModEq{T}(X::ModEq{T}) where {T} = new{T}(X.p, X.f, X.g, X.h, X.k, X.L)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
ModEq(X::AbstractVector{T}) where {T} = ModEq{T}(X[1], X[2], X[3], X[4], X[5], X[6])
function ModEq(p::P, f::F, g::G, h::H, k::K, L::LT) where {P,F,G,H,K,LT}
    return ModEq{promote_type(P, F, G, H, K, LT)}(p, f, g, h, k, L)
end
# More specific than AbstractVector to avoid ambiguity
ModEq(g::StaticVector{N,T}) where {N,T} = ModEq{T}(g[1], g[2], g[3], g[4], g[5], g[6])
ModEq{T}(g::StaticVector) where {T} = ModEq{T}(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::ModEq{T}) where {T<:Number} = SVector{6,T}(g.p, g.f, g.g, g.h, g.k, g.L)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{M}; T::DataType=Float64) where {M<:ModEq}
    return M{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::ModEq{T}, i::Int) where {T<:Number}
    if i == 1
        return p.p
    elseif i == 2
        return p.f
    elseif i == 3
        return p.g
    elseif i == 4
        return p.h
    elseif i == 5
        return p.k
    elseif i == 6
        return p.L
    else
        throw(BoundsError(p, i))
    end
end

export ModEqN
"""
    ModEqN{T} <: AstroCoord

Modified Equinoctial Orbital Elements with Mean Motion. 6D parameterziation of the orbit.

# Fields
- `η` - mean motion (rad/s)
- `f` - eccentricity projection onto longitude of perigee
- `g` - eccentricity projection onto ⟂ longitude of perigee
- `h` - projection of half inclination onto RAAN
- `k` - projection of half inclination onto ⟂ RAAN
- `L` - true longitude

# Constructors
- `ModEqN(η, f, g, h, k, L)`
- `ModEqN(X::AbstractArray)`
- `ModEqN(X::AstroCoord, μ::Number)`

# Examples
```julia
# Create ModEqN coordinates with mean motion
μ = 3.986004418e14  # Earth's gravitational parameter
a = 7000e3          # Semi-major axis
η = √(μ / a^3)      # Mean motion
f = 0.1             # Eccentricity projection

modeqn = ModEqN(η, f, 0.0, 0.0, 0.0, π/3)
modeqn.η  # 0.001078 (mean motion)
modeqn.f  # 0.1 (f component)
```

# See Also
- [`ModEq`](@ref): Modified Equinoctial coordinates with semi-parameter instead of mean motion
"""
struct ModEqN{T} <: AstroCoord{6,T}
    η::T
    f::T
    g::T
    h::T
    k::T
    L::T
    @inline ModEqN{T}(η, f, g, h, k, L) where {T} = new{T}(η, f, g, h, k, L)
    @inline ModEqN{T}(X::ModEqN{T}) where {T} = new{T}(X.η, X.f, X.g, X.h, X.k, X.L)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
ModEqN(X::AbstractVector{T}) where {T} = ModEqN{T}(X[1], X[2], X[3], X[4], X[5], X[6])
function ModEqN(η::ETA, f::F, g::G, h::H, k::K, L::LT) where {ETA,F,G,H,K,LT}
    return ModEqN{promote_type(ETA, F, G, H, K, LT)}(η, f, g, h, k, L)
end
# More specific than AbstractVector to avoid ambiguity
ModEqN(g::StaticVector{N,T}) where {N,T} = ModEqN{T}(g[1], g[2], g[3], g[4], g[5], g[6])
ModEqN{T}(g::StaticVector) where {T} = ModEqN{T}(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::ModEqN{T}) where {T<:Number} = SVector{6,T}(g.η, g.f, g.g, g.h, g.k, g.L)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{M}; T::DataType=Float64) where {M<:ModEqN}
    return M{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::ModEqN{T}, i::Int) where {T<:Number}
    if i == 1
        return p.η
    elseif i == 2
        return p.f
    elseif i == 3
        return p.g
    elseif i == 4
        return p.h
    elseif i == 5
        return p.k
    elseif i == 6
        return p.L
    else
        throw(BoundsError(p, i))
    end
end

export J2EqOE
"""
    J2EqOE{T} <: AstroCoord

Modified Equinoctial Orbital Elements. 6D parameterziation of the orbit.
n - mean motion
h - eccetricity projection onto longitude of perigee
k - eccetricity projection onto ⟂ longitude of perigee
p - projection of half inclination onto RAAN
q - projection of half inclination onto ⟂ RAAN
L - true longitude

Constructors
J2EqOE(n, h, k, p, q, L)
J2EqOE(X::AbstractArray)
J2EqOE(X::AstroCoord, μ::Number)

"""
struct J2EqOE{T} <: AstroCoord{6,T}
    n::T
    h::T
    k::T
    p::T
    q::T
    L::T
    @inline J2EqOE{T}(n, h, k, p, q, L) where {T} = new{T}(n, h, k, p, q, L)
    @inline J2EqOE{T}(X::J2EqOE{T}) where {T} = new{T}(X.n, X.h, X.k, X.p, X.q, X.L)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
J2EqOE(X::AbstractVector{T}) where {T} = J2EqOE{T}(X[1], X[2], X[3], X[4], X[5], X[6])
function J2EqOE(n::N, h::H, k::K, p::P, q::Q, L::LT) where {N,H,K,P,Q,LT}
    return J2EqOE{promote_type(N, H, K, P, Q, LT)}(n, h, k, p, q, L)
end
# More specific than AbstractVector to avoid ambiguity
J2EqOE(g::StaticVector{N,T}) where {N,T} = J2EqOE{T}(g[1], g[2], g[3], g[4], g[5], g[6])
J2EqOE{T}(g::StaticVector) where {T} = J2EqOE{T}(g[1], g[2], g[3], g[4], g[5], g[6])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::J2EqOE{T}) where {T<:Number} = SVector{6,T}(g.n, g.h, g.k, g.p, g.q, g.L)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
function Base.one(::Type{J}; T::DataType=Float64) where {J<:J2EqOE}
    return J{T}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function Base.getindex(p::J2EqOE{T}, i::Int) where {T<:Number}
    if i == 1
        return p.n
    elseif i == 2
        return p.h
    elseif i == 3
        return p.k
    elseif i == 4
        return p.p
    elseif i == 5
        return p.q
    elseif i == 6
        return p.L
    else
        throw(BoundsError(p, i))
    end
end
