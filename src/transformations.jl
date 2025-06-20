"""
    abstract type AstrodynamicsTransformation <: Transformation

An abstract type representing a Transformation Between Astrodynamics Coordinates
"""
abstract type AstrodynamicsTransformation <: Transformation end

"""
    abstract type AstroCoordTransformation <: AstrodynamicsTransformation

An abstract type representing a Transformation Between Astrodynamics Coordinates with no Time Regularization
"""
abstract type AstroCoordTransformation <: AstrodynamicsTransformation end

# ~~~~~~~~~~~~~~~ Transformation Definitions ~~~~~~~~~~~~~~~ #
macro define_transformation_pair(
    BaseFrame, TargetFrame, base_to_target_backend, target_to_base_backend
)
    BaseFrameType = esc(BaseFrame)
    TargetFrameType = esc(TargetFrame)
    base_to_target_transform = esc(Symbol(BaseFrame, :To, TargetFrame, :Transform))
    base_to_target = esc(Symbol(BaseFrame, :To, TargetFrame))
    target_to_base_transform = esc(Symbol(TargetFrame, :To, BaseFrame, :Transform))
    target_to_base = esc(Symbol(TargetFrame, :To, BaseFrame))

    quote
        struct $base_to_target_transform <: AstroCoordTransformation end
        @inline function (::$base_to_target_transform)(
            x::$BaseFrameType{T}, μ::V
        ) where {T<:Number,V<:Number}
            RT = promote_type(T, V)
            return $TargetFrameType{RT}($(esc(base_to_target_backend))(params(x), μ))
        end
        const $base_to_target = $base_to_target_transform()

        struct $target_to_base_transform <: AstroCoordTransformation end
        @inline function (::$target_to_base_transform)(
            x::$TargetFrameType{T}, μ::V
        ) where {T<:Number,V<:Number}
            RT = promote_type(T, V)
            return $BaseFrameType{RT}($(esc(target_to_base_backend))(params(x), μ))
        end
        const $target_to_base = $target_to_base_transform()

        Base.inv(::$base_to_target_transform) = $target_to_base_transform()
        Base.inv(::$target_to_base_transform) = $base_to_target_transform()
    end
end

# ~~~~~~~~~~~~~~~ Direct Transformations ~~~~~~~~~~~~~~~ #
@define_transformation_pair Cartesian Keplerian cart2koe koe2cart
@define_transformation_pair Cartesian Milankovich cart2Mil Mil2cart
@define_transformation_pair Cartesian Cylindrical cart2cylind cylind2cart
@define_transformation_pair Cartesian Spherical cart2sphere sphere2cart
@define_transformation_pair Cartesian Delaunay cart2delaunay delaunay2cart
@define_transformation_pair Cartesian J2EqOE cart2J2EqOE J2EqOE2cart
@define_transformation_pair Keplerian USM7 koe2USM7 USM72koe
@define_transformation_pair Keplerian ModEq koe2ModEq ModEq2koe
@define_transformation_pair USM7 USM6 USM72USM6 USM62USM7
@define_transformation_pair USM7 USMEM USM72USMEM USMEM2USM7

# ~~~~~~~~~~~~~~~ Aliases for Readability ~~~~~~~~~~~~~~~ #
const KeplerianToModifiedEquinoctial = KeplerianToModEq
const ModifiedEquinoctialToKeplerian = ModEqToKeplerian

# ~~~~~~~~~~~~~~~ All Composed Transformations ~~~~~~~~~~~~~~~ #
const COORD_TYPES = (
    Cartesian,
    Keplerian,
    USM7,
    USM6,
    USMEM,
    Milankovich,
    ModEq,
    Cylindrical,
    Spherical,
    Delaunay,
    J2EqOE,
)
const COORD_NAMES = Dict(
    Cartesian => :Cartesian,
    Keplerian => :Keplerian,
    USM7 => :USM7,
    USM6 => :USM6,
    USMEM => :USMEM,
    Milankovich => :Milankovich,
    ModEq => :ModifiedEquinoctial,
    Cylindrical => :Cylindrical,
    Spherical => :Spherical,
    Delaunay => :Delaunay,
    J2EqOE => :J2EqOE,
)

# Build a graph of transformations to find paths
TRANSFORM_GRAPH = Dict{Symbol,Vector{Symbol}}()
for T in COORD_TYPES
    T_name = COORD_NAMES[T]
    TRANSFORM_GRAPH[T_name] = []
end

function add_transform_edge(T1, T2)
    T1_name = COORD_NAMES[T1]
    T2_name = COORD_NAMES[T2]
    push!(TRANSFORM_GRAPH[T1_name], T2_name)
    push!(TRANSFORM_GRAPH[T2_name], T1_name)
end

add_transform_edge(Cartesian, Keplerian)
add_transform_edge(Cartesian, Milankovich)
add_transform_edge(Cartesian, Cylindrical)
add_transform_edge(Cartesian, Spherical)
add_transform_edge(Cartesian, Delaunay)
add_transform_edge(Cartesian, J2EqOE)
add_transform_edge(Keplerian, USM7)
add_transform_edge(Keplerian, ModEq)
add_transform_edge(USM7, USM6)
add_transform_edge(USM7, USMEM)

# Use BFS to find transformation paths between any two coordinate systems
function find_transform_path(start_node, end_node)
    queue = [[start_node]]
    visited = Set([start_node])

    while !isempty(queue)
        path = popfirst!(queue)
        node = path[end]

        if node == end_node
            return path
        end

        for neighbor in TRANSFORM_GRAPH[node]
            if !(neighbor in visited)
                push!(visited, neighbor)
                new_path = copy(path)
                push!(new_path, neighbor)
                push!(queue, new_path)
            end
        end
    end
    return nothing # Should not happen in a connected graph
end

# Define all transformations based on paths
for T1 in COORD_TYPES
    for T2 in COORD_TYPES
        T1 == T2 && continue

        T1_name = COORD_NAMES[T1]
        T2_name = COORD_NAMES[T2]

        transform_name = Symbol(T1_name, :To, T2_name)
        isdefined(@__MODULE__, transform_name) && continue

        path = find_transform_path(T1_name, T2_name)

        if path !== nothing && length(path) > 1
            # Chain the transformations together
            transforms_to_compose = []
            for i in 1:(length(path) - 1)
                from = path[i]
                to = path[i + 1]
                t_name = Symbol(from, :To, to)
                push!(transforms_to_compose, getfield(@__MODULE__, t_name))
            end

            # Compose from right to left
            composed_transform = transforms_to_compose[end]
            for i in (length(transforms_to_compose) - 1):-1:1
                composed_transform = composed_transform ∘ transforms_to_compose[i]
            end

            @eval const $transform_name = $composed_transform

            # Define inverse automatically
            inv_transform_name = Symbol(T2_name, :To, T1_name)
            isdefined(@__MODULE__, inv_transform_name) && continue
            @eval const $inv_transform_name = inv($transform_name)
        end
    end
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~ Additional Constructors ~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
for ToCoord in COORD_TYPES
    ToCoordName = ToCoord isa UnionAll ? ToCoord.body.name.name : ToCoord.name.name
    for FromCoord in COORD_TYPES
        if ToCoord == FromCoord
            @eval ($ToCoordName)(X::$ToCoord{T}, μ::Number) where {T<:Number} = X
        else
            ToCoordTransformName = COORD_NAMES[ToCoord]
            FromCoordTransformName = COORD_NAMES[FromCoord]
            Transform = Symbol(FromCoordTransformName, :To, ToCoordTransformName)
            @eval ($ToCoordName)(X::$FromCoord{T}, μ::Number) where {T<:Number} = $(getfield(
                @__MODULE__, Transform
            ))(
                X, μ
            )
        end
    end
end
