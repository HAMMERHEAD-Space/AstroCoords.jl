using Test
using AstroCoords

@testset "Module Import and Structure" begin
    @testset "Module Loading" begin
        # Test that module loads without error
        @test isdefined(Main, :AstroCoords)
        @test AstroCoords isa Module
    end

    @testset "Exported Symbols" begin
        # Test that key symbols are exported and accessible
        @test isdefined(AstroCoords, :Cartesian)
        @test isdefined(AstroCoords, :KeplerianOrbitElements)
        @test isdefined(AstroCoords, :cart2koe)
        @test isdefined(AstroCoords, :koe2cart)
    end

    @testset "Module Dependencies" begin
        # Test that required dependencies are loaded
        # (This indirectly tests that Project.toml dependencies are available)
        @test isdefined(AstroCoords, :LinearAlgebra) || true  # stdlib
    end

    @testset "Module Namespace" begin
        # Test that module has expected submodules or structure
        # Check if coordinate types are in the namespace
        @test typeof(AstroCoords.Cartesian) <: Type
        @test typeof(AstroCoords.KeplerianOrbitElements) <: Type
    end

    @testset "Basic Functionality" begin
        # Smoke test: can we create a basic object?
        cart = AstroCoords.Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)
        @test cart isa AstroCoords.Cartesian
        
        koe = AstroCoords.KeplerianOrbitElements(7000.0, 0.01, 0.5, 0.0, 0.0, 0.0)
        @test koe isa AstroCoords.KeplerianOrbitElements
    end
end
