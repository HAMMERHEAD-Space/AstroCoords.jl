module AstroCoords

using LinearAlgebra
using StaticArrays

include("utils.jl")
include("./anomalies.jl")
include("./attitude_changes.jl")
include("./coordinate_changes/coordinate_changes.jl")
include("./coordinate_changes/J2EqOE.jl")
include("./coordinate_changes/EDromo.jl")
include("./core_types.jl")

include("./coordinate_sets/cartesian.jl")
include("./coordinate_sets/delaunay.jl")
include("./coordinate_sets/keplerian.jl")
include("./coordinate_sets/milankovich.jl")
include("./coordinate_sets/modEq.jl")
include("./coordinate_sets/spherical.jl")
include("./coordinate_sets/usm.jl")
include("./coordinate_sets/edromo.jl")

include("./transformations.jl")

include("./quantities.jl")

end
