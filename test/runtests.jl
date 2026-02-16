using AstroCoords
using LinearAlgebra
using StaticArrays
using Test

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
    Poincare,
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
    include("test_GEqOE.jl")
end

# ── Gated differentiability tests ────────────────────────────────────────────
const _DIFF_ENV = get(ENV, "ASTROCOORDS_TEST_DIFF", "false")
if _DIFF_ENV ∉ ("false", "")
    using DifferentiationInterface
    using FiniteDiff

    _run_all = _DIFF_ENV ∈ ("true", "all")
    _requested = _run_all ? Set{String}() : Set(strip.(split(_DIFF_ENV, ",")))
    _need(name) = _run_all || name ∈ _requested

    _backend_list = Tuple{String,Any}[]

    if _need("ForwardDiff")
        using ForwardDiff
        push!(_backend_list, ("ForwardDiff", AutoForwardDiff()))
    end
    if _need("Enzyme")
        using Enzyme
        push!(_backend_list, ("Enzyme", AutoEnzyme()))
    end
    if _need("Mooncake")
        using Mooncake
        push!(_backend_list, ("Mooncake", AutoMooncake(; config=nothing)))
    end
    if _need("PolyesterForwardDiff")
        using PolyesterForwardDiff
        push!(_backend_list, ("PolyesterForwardDiff", AutoPolyesterForwardDiff()))
    end

    const _TEST_ZYGOTE = _need("Zygote")
    if _TEST_ZYGOTE
        using Zygote
    end

    isempty(_backend_list) &&
        !_TEST_ZYGOTE &&
        error("No matching backend for ASTROCOORDS_TEST_DIFF=\"$_DIFF_ENV\"")

    const _BACKENDS = Tuple(_backend_list)
    @info "Running differentiability tests with: $(join([b[1] for b in _BACKENDS], ", "))$(_TEST_ZYGOTE ? ", Zygote" : "")"

    @testset "Differentiation" begin
        include("test_differentiability.jl")
    end
else
    @info "Skipping differentiability tests (set ASTROCOORDS_TEST_DIFF to enable)"
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
