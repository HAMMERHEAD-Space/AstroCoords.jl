export Keplerian
"""
    Keplerian{T} <: AstroCoord

Keplerian Orbital Elements. 6D parameterziation of the orbit.
a - semi-major axis
e - eccetricity 
i - inclination 
Ω - Right Ascension of Ascending Node 
ω - Argument of Perigee
f - True Anomaly

Constructors
Keplerian(a, e, i, Ω, ω, f)
Keplerian(X::AbstractArray)
Keplerian(X::AstroCoord, μ::Number)

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

#TODO: Return other anomalies
