module AstroCoords

using LinearAlgebra
using StaticArrays

# attitude
export EP2MRP, MRP2EP
export DCM2EP, EP2DCM
export DCM2MRP, MRP2DCM

# coordinate sets
export CoordinateSet
export Cartesian
export Keplerian
export ModEq
export Delaunay
export Milankovich
export Spherical
export Cylindrical
export USM6, USM7, USMEM
export J2EqOE
export EDromo, KustaanheimoStiefel, StiefelScheifele

export params
export orbitalNRG
export angularMomentumQuantity

# regularized coordinates
export RegularizedCoordinateConfig
export PhysicalTime, ConstantTime, LinearTime
export compute_characteristic_scales
export compute_initial_phi

include("angle_utils.jl")
include("anomalies.jl")
include("quantities.jl")
include("attitude/attitude.jl")
include("coordinate_sets.jl")
include("regularized/regularized.jl")

const COORD_TYPES = (
    Cartesian,
    Keplerian,
    ModEq,
    Delaunay,
    Milankovich,
    Spherical,
    Cylindrical,
    USM6,
    USM7,
    USMEM,
    J2EqOE,
    EDromo,
    KustaanheimoStiefel,
    StiefelScheifele,
)

end
