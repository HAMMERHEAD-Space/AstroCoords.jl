using AstroCoords
using LinearAlgebra
using StaticArrays
using Test

using DifferentiationInterface
using Diffractor
using Enzyme
using FiniteDiff
using ForwardDiff
using Mooncake
using PolyesterForwardDiff
using Zygote

using Aqua
using JET
using AllocCheck

const _COORDINATE_SETS = [
    Cartesian,
    Delaunay,
    Keplerian,
    Milankovich,
    ModEq,
    ModEqN,
    Cylindrical,
    Spherical,
    USM7,
    USM6,
    USMEM,
    J2EqOE,
    EDromo,
    KustaanheimoStiefel,
    StiefelScheifele,
]

@testset "AstroCoords.jl" begin
    # Core Types
    include("test_core_types.jl")
    include("test_regularized_config.jl")
    include("test_transformations.jl")

    # Coordinate Constructors
    include("coordinate_sets/test_cartesian_constructors.jl")
    include("coordinate_sets/test_delaunay.jl")
    include("coordinate_sets/test_keplerian.jl")
    include("coordinate_sets/test_ks_set.jl")
    include("coordinate_sets/test_milankovich.jl")
    include("coordinate_sets/test_modeq.jl")
    include("coordinate_sets/test_spherical.jl")
    include("coordinate_sets/test_ss_set.jl")
    include("coordinate_sets/test_usm.jl")

    # Utility Functions
    include("test_anomalies.jl")
    include("test_quantities.jl")
    include("test_utils.jl")

    # Coordinate Changes
    include("test_attitude_changes.jl")
    include("test_coordinate_changes.jl")
    include("test_J2EqOE.jl")
    include("test_EDromo.jl")
    include("test_KS.jl")
    include("test_StiefelScheifele.jl")
end

@testset "Differentiation" begin
    include("test_differentiability.jl")
end

@testset "Code Performance" begin
    include("test_JET.jl")
    include("test_allocs.jl")
end

@testset "Aqua.jl" begin
    Aqua.test_all(
        AstroCoords; ambiguities=(recursive = false), deps_compat=(check_extras = false)
    )
end
