using AstroCoords
using Test
using AllocCheck
using DifferentiationInterface
using FiniteDiff
using ForwardDiff
using Enzyme
using Diffractor
using Mooncake
using PolyesterForwardDiff
using Zygote
using JET
using LinearAlgebra
using StaticArrays

const _COORDINATE_SETS = (
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

@testset "AstroCoords.jl" begin
    include("test_anomalies.jl")
    include("test_quantities.jl")
    include("test_coordinate_changes.jl")
    include("test_J2EqOE.jl")
    include("test_EDromo.jl")
    include("test_KS.jl")
    include("test_StiefelScheifele.jl")
    include("test_differentiability.jl")
    include("test_allocs.jl")
    include("test_JET.jl")
    include("test_coverage.jl")
end
