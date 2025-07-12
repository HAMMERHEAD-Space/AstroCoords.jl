export Coordinate
"""
    abstract type AstroCoord{N, T} <: StaticMatrix{N, 1, T}

An abstract type representing a N-Dimensional Coordinate Set
"""
abstract type Coordinate{N,T} <: StaticVector{N,T} end

Base.@pure StaticArrays.Size(::Type{Coordinate{N}}) where {N} = Size(N)
Base.@pure StaticArrays.Size(::Type{Coordinate{N,T}}) where {N,T} = Size(N)
Base.@pure StaticArrays.Size(::Type{A}) where {A<:Coordinate} = Size(supertype(A))

# Generate zero-matrix with SVector
# Note that zeros(AstroCoord,dims...) is not Array{<:AstroCoord} but Array{<:StaticArray{N}}
Base.zero(::Coordinate{N,T}) where {N,T} = @SVector zeros(T, N)
Base.zero(::Type{Coordinate}) = error("The dimension of the Coordinate is not specified.")
Base.zero(::Type{<:Coordinate{N}}) where {N} = @SVector zeros(N)
Base.zero(::Type{<:Coordinate{N,T}}) where {N,T} = @SVector zeros(T, N)
Base.zeros(::Type{A}) where {A<:Coordinate} = zeros(A) # avoid StaticArray constructor
function Base.zeros(::Type{A}, dims::Base.DimOrInd...) where {A<:Coordinate}
    return zeros(typeof(zero(A)), dims...)
end
function Base.zeros(::Type{A}, dims::NTuple{N,Integer}) where {A<:Coordinate,N}
    return zeros(typeof(zero(A)), dims)
end
Base.zeros(::Type{A}, dims::Tuple{}) where {A<:Coordinate} = zeros(typeof(zero(A)), dims) # avoid ambiguity

# `convert` goes through the constructors, similar to e.g. `Number`
Base.convert(::Type{A}, coord::A) where {N,A<:Coordinate{N}} = coord
Base.convert(::Type{A}, coord::Coordinate{N}) where {N,A<:Coordinate{N}} = A(coord)

abstract type AstroCoord{N,T} <: Coordinate{N,T} end
export params
abstract type AttitudeCoord{N,T} <: Coordinate{N,T} end

"""
The `Transformation` supertype defines a simple interface for performing
transformations. Subtypes should be able to apply a coordinate system
transformation on the correct data types by overloading the call method, and
usually would have the corresponding inverse transformation defined by `Base.inv()`.
Efficient compositions can optionally be defined by `compose()` (equivalently `∘`).
"""
abstract type Transformation end

"""
The `IdentityTransformation` is a singleton `Transformation` that returns the
input unchanged, similar to `identity`.
"""
struct IdentityTransformation <: Transformation end

@inline (::IdentityTransformation)(x) = x

struct ComposedTransformation{T1<:Transformation,T2<:Transformation} <: Transformation
    t1::T1
    t2::T2
end

@inline function (trans::ComposedTransformation)(x::AstroCoord, μ::Number; kwargs...)
    return trans.t1(trans.t2(x, μ; kwargs...), μ; kwargs...)
end

function compose(trans1::T, trans2::V) where {T<:Transformation} where {V<:Transformation}
    return ComposedTransformation(trans1, trans2)
end

function Base.:∘(trans1::T, trans2::V) where {T<:Transformation} where {V<:Transformation}
    return ComposedTransformation(trans1, trans2)
end

compose(trans::IdentityTransformation, ::IdentityTransformation) = trans
compose(::IdentityTransformation, trans::Transformation) = trans
compose(trans::Transformation, ::IdentityTransformation) = trans

"""
    inv(trans::Transformation)
Returns the inverse (or reverse) of the transformation `trans`
"""
function Base.inv(trans::Transformation)
    return error("Inverse transformation for $(typeof(trans)) has not been defined.")
end

Base.inv(trans::ComposedTransformation) = inv(trans.t2) ∘ inv(trans.t1)
Base.inv(trans::IdentityTransformation) = trans
