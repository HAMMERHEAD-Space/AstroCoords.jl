module AstroCoords

using LinearAlgebra
using StaticArrays

include("./core_types.jl")
include("utils.jl")
include("./anomalies.jl")
include("./attitude_changes.jl")

include("./coordinate_changes/coordinate_changes.jl")
include("./coordinate_changes/J2EqOE.jl")
include("./coordinate_changes/EDromo.jl")
include("./coordinate_changes/Kustaanheimo-Stiefel.jl")
include("./coordinate_changes/Stiefel-Scheifele.jl")

include("./coordinate_sets/cartesian.jl")
include("./coordinate_sets/delaunay.jl")
include("./coordinate_sets/keplerian.jl")
include("./coordinate_sets/milankovich.jl")
include("./coordinate_sets/modEq.jl")
include("./coordinate_sets/spherical.jl")
include("./coordinate_sets/usm.jl")
include("./coordinate_sets/edromo.jl")
include("./coordinate_sets/kustaanheimo-stiefel.jl")
include("./coordinate_sets/stiefel-scheifele.jl")

include("./transformations.jl")

const COORDINATE_SET_ALIASES = Dict(
    "Cartesian" => Cartesian,
    "Cylindrical" => Cylindrical,
    "Delaunay" => Delaunay,
    "EDromo" => EDromo,
    "J2EqOE" => J2EqOE,
    "Keplerian" => Keplerian,
    "KustaanheimoStiefel" => KustaanheimoStiefel,
    "StiefelScheifele" => StiefelScheifele,
    "Milankovich" => Milankovich,
    "ModifiedEquinoctial" => ModEq,
    "Spherical" => Spherical,
    "USM6" => USM6,
    "USM7" => USM7,
    "USMEM" => USMEM,
)

include("./quantities.jl")

end
