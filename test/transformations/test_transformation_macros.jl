using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "Transformation Macro System" begin
    μ = 398600.4418  # Earth's gravitational parameter (km^3/s^2)

    @testset "Transformation Struct Creation" begin
        # Tests lines 16-46 - macro expansion creates proper structs
        @test isdefined(AstroCoords, :CartesianToKeplerianTransform)
        @test isdefined(AstroCoords, :KeplerianToCartesianTransform)
        
        # Check singleton instances
        @test isdefined(AstroCoords, :CartesianToKeplerian)
        @test isdefined(AstroCoords, :KeplerianToCartesian)
        
        # Verify types
        @test AstroCoords.CartesianToKeplerian isa AstroCoords.CartesianToKeplerianTransform
        @test AstroCoords.KeplerianToCartesian isa AstroCoords.KeplerianToCartesianTransform
    end

    @testset "Forward Transformation: Cartesian → Keplerian" begin
        # Tests lines 71-191 - actual transformation dispatch
        r = 7000.0
        v = sqrt(μ / r)
        cart = Cartesian(r, 0.0, 0.0, 0.0, v, 0.0)
        
        kep = AstroCoords.CartesianToKeplerian(cart, μ)
        
        @test kep isa Keplerian
        @test kep.a ≈ r atol=1e-10
        @test kep.e < 1e-14  # Circular orbit
    end

    @testset "Inverse Transformation: Keplerian → Cartesian" begin
        # Tests lines 164-196 - inverse transformations work
        a = 8000.0
        e = 0.15
        i = 0.5
        Ω = 1.0
        ω = 0.5
        f = 0.3
        
        kep = Keplerian(a, e, i, Ω, ω, f)
        cart = AstroCoords.KeplerianToCartesian(kep, μ)
        
        @test cart isa Cartesian
        @test norm(SVector{3}(cart.x, cart.y, cart.z)) ≈ a*(1-e^2)/(1+e*cos(f)) atol=1e-10
    end

    @testset "Inverse Function Behavior" begin
        # Tests lines 164-196 - inv() returns correct inverse transform
        forward = AstroCoords.CartesianToKeplerianTransform()
        backward = inv(forward)
        
        @test backward isa AstroCoords.KeplerianToCartesianTransform
        
        # Double inverse should return original
        forward_again = inv(backward)
        @test forward_again isa AstroCoords.CartesianToKeplerianTransform
    end

    @testset "Caching/Singleton Behavior" begin
        # Tests lines 235-244 - verify singleton pattern
        t1 = AstroCoords.CartesianToKeplerian
        t2 = AstroCoords.CartesianToKeplerian
        
        # Should be the exact same object (singleton)
        @test t1 === t2
    end

    @testset "Transformation Composition: A→B→C" begin
        # Tests lines 262-323 - composing transformations
        r = 7500.0
        v = sqrt(μ / r) * 0.95
        cart_orig = Cartesian(r, 0.0, 0.0, 0.0, v, v*0.1)
        
        # Cart → Kep → ModEq → Kep → Cart
        kep = AstroCoords.CartesianToKeplerian(cart_orig, μ)
        modEq = AstroCoords.KeplerianToModEq(kep, μ)
        kep_restored = AstroCoords.ModEqToKeplerian(modEq, μ)
        cart_restored = AstroCoords.KeplerianToCartesian(kep_restored, μ)
        
        @test cart_restored.x ≈ cart_orig.x atol=1e-10 rtol=1e-12
        @test cart_restored.y ≈ cart_orig.y atol=1e-10 rtol=1e-12
        @test cart_restored.z ≈ cart_orig.z atol=1e-10 rtol=1e-12
        @test cart_restored.ẋ ≈ cart_orig.ẋ atol=1e-10 rtol=1e-12
        @test cart_restored.ẏ ≈ cart_orig.ẏ atol=1e-10 rtol=1e-12
        @test cart_restored.ż ≈ cart_orig.ż atol=1e-10 rtol=1e-12
    end

    @testset "Multiple Transformation Pairs" begin
        # Verify multiple defined transformations
        @test isdefined(AstroCoords, :CartesianToMilankovichTransform)
        @test isdefined(AstroCoords, :CartesianToCylindricalTransform)
        @test isdefined(AstroCoords, :CartesianToSphericalTransform)
        @test isdefined(AstroCoords, :CartesianToDelaunayTransform)
        @test isdefined(AstroCoords, :KeplerianToUSM7Transform)
        @test isdefined(AstroCoords, :KeplerianToModEqTransform)
        @test isdefined(AstroCoords, :USM7ToUSM6Transform)
    end

    @testset "Transformation with Different Types" begin
        # Test type promotion in transformations
        cart_f32 = Cartesian(Float32(7000), 0.0f0, 0.0f0, 0.0f0, Float32(sqrt(398600.4418/7000)), 0.0f0)
        kep = AstroCoords.CartesianToKeplerian(cart_f32, Float64(μ))
        
        @test eltype(params(kep)) == Float64  # Should promote to Float64
    end

    @testset "Alias Definitions" begin
        # Test that aliases work (lines 235-244)
        @test AstroCoords.KeplerianToModifiedEquinoctial === AstroCoords.KeplerianToModEq
        @test AstroCoords.ModifiedEquinoctialToKeplerian === AstroCoords.ModEqToKeplerian
    end

    @testset "Coordinate Type Interface" begin
        # Verify that all coordinate types work with transformations
        r = 7000.0
        v = sqrt(μ / r)
        cart = Cartesian(r, 0.0, 0.0, 0.0, v, 0.0)
        
        # Test params() extraction
        @test length(params(cart)) == 6
        @test params(cart)[1] == r
        
        # Test construction from vector
        kep = AstroCoords.CartesianToKeplerian(cart, μ)
        kep_params = params(kep)
        @test length(kep_params) == 6
    end

    @testset "USM Transformation Chain" begin
        # Test USM7 ↔ USM6 and USM7 ↔ USMEM transformations
        a = 8000.0
        e = 0.2
        i = 0.5
        Ω = 1.0
        ω = 0.5
        f = 0.3
        
        kep = Keplerian(a, e, i, Ω, ω, f)
        usm7 = AstroCoords.KeplerianToUSM7(kep, μ)
        
        @test usm7 isa USM7
        
        # Test USM7 → USM6
        usm6 = AstroCoords.USM7ToUSM6(usm7, μ)
        @test usm6 isa USM6
        
        # Roundtrip: USM7 → USM6 → USM7
        usm7_restored = AstroCoords.USM6ToUSM7(usm6, μ)
        @test usm7_restored isa USM7
    end
end
