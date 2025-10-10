module AstroCoords

using LinearAlgebra
using StaticArrays

export Cartesian,
    Cylindrical,
    Spherical,
    Keplerian,
    ModEq,
    Milankovich,
    Delaunay,
    USM6,
    USM7,
    USMEM,
    J2EqOE,
    EDromo,
    KustaanheimoStiefel,
    StiefelScheifele

export PhysicalTime, ConstantTime, LinearTime
export RegularizedCoordinateConfig

# Include base types first
include("coordinate_sets/coordinate_set.jl")

# Then include all other coordinate set files
for file in readdir(joinpath(@__DIR__, "coordinate_sets"); join=true)
    # Skip the base file since we already included it
    if !endswith(file, "coordinate_set.jl")
        include(file)
    end
end

include("quantities.jl")
include("attitude.jl")
include("anomalies.jl")
include("coordinate_changes.jl")
include("utils.jl")

const COORD_TYPES = (
    Cartesian,
    Cylindrical,
    Delaunay,
    J2EqOE,
    Keplerian,
    Milankovich,
    ModEq,
    Spherical,
    USM6,
    USM7,
    USMEM,
)

end
