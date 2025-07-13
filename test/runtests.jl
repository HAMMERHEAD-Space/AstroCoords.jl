using AstroCoords
using Aqua
using Test

using DifferentiationInterface
using Diffractor
using Enzyme
using FiniteDiff
using ForwardDiff
using Mooncake
using PolyesterForwardDiff
using Zygote

using JET
using AllocCheck

const _COORDINATE_SETS = [
    Cartesian,
    Delaunay,
    Keplerian,
    Milankovich,
    ModEq,
    Cylindrical,
    Spherical,
    USM7,
    USM6,
    USMEM,
    J2EqOE,
    EDromo,
]

@testset "AstroCoords.jl" begin
    include("test_J2EqOE.jl")
    include("test_EDromo.jl")
    include("test_coordinate_changes.jl")
    include("test_anomalies.jl")
    include("test_quantities.jl")
end

@testset "Differentiation" begin
    include("test_differentiability.jl")
end

@testset "Code Performance" begin
    include("test_JET.jl")
    include("test_allocs.jl")
end

@testset "Aqua.jl" begin
    Aqua.test_all(AstroCoords; ambiguities=(recursive = false))
end
