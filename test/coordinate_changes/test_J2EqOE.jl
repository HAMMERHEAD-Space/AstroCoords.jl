using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "J2EqOE Transformations" begin
    μ = 398600.4418  # Earth's gravitational parameter (km^3/s^2)

    @testset "koe2IOE: Standard Elliptical Orbit" begin
        # Standard elliptical orbit
        a = 7000.0  # km
        e = 0.01
        i = 0.5  # rad (~28.6 degrees)
        Ω = 1.0  # rad
        ω = 0.5  # rad
        f = 0.3  # rad (true anomaly)
        
        koe = SVector{6}(a, e, i, Ω, ω, f)
        IOE = AstroCoords.koe2IOE(koe, μ)
        
        # IOE should have 6 elements
        @test length(IOE) == 6
        @test IOE[1] ≈ a atol=1e-10  # I1 = a
        @test norm(SVector{2}(IOE[2], IOE[3])) ≈ e atol=1e-10  # e = sqrt(I2^2 + I3^2)
        @test 2.0 * asin(sqrt(IOE[4]^2 + IOE[5]^2)) ≈ i atol=1e-10  # i = 2*asin(sqrt(I4^2+I5^2))
    end

    @testset "IOE2koe: Roundtrip Conversion" begin
        # Test roundtrip: koe → IOE → koe
        a = 8000.0
        e = 0.15
        i = 1.2
        Ω = 2.5
        ω = 1.8
        f = 0.7
        
        koe_orig = SVector{6}(a, e, i, Ω, ω, f)
        IOE = AstroCoords.koe2IOE(koe_orig, μ)
        koe_restored = AstroCoords.IOE2koe(IOE, μ)
        
        @test koe_restored[1] ≈ a atol=1e-10 rtol=1e-12
        @test koe_restored[2] ≈ e atol=1e-10 rtol=1e-12
        @test koe_restored[3] ≈ i atol=1e-10 rtol=1e-12
        @test koe_restored[4] ≈ Ω atol=1e-10 rtol=1e-12
        @test koe_restored[5] ≈ ω atol=1e-10 rtol=1e-12
        @test koe_restored[6] ≈ f atol=1e-10 rtol=1e-12
    end

    @testset "koe2IOE: Circular Orbit (e≈0)" begin
        a = 7000.0
        e = 1e-10  # Nearly circular
        i = 0.5
        Ω = 1.0
        ω = 0.5
        f = 0.3
        
        koe = SVector{6}(a, e, i, Ω, ω, f)
        IOE = AstroCoords.koe2IOE(koe, μ)
        
        @test IOE[1] ≈ a atol=1e-10
        @test norm(SVector{2}(IOE[2], IOE[3])) ≈ e atol=1e-14  # Very small eccentricity
    end

    @testset "koe2IOE: Near-Parabolic Orbit (e≈1)" begin
        a = 7000.0
        e = 0.9999  # Nearly parabolic
        i = 0.5
        Ω = 1.0
        ω = 0.5
        f = 0.1  # Small true anomaly to avoid issues
        
        koe = SVector{6}(a, e, i, Ω, ω, f)
        IOE = AstroCoords.koe2IOE(koe, μ)
        
        @test IOE[1] ≈ a atol=1e-10
        @test norm(SVector{2}(IOE[2], IOE[3])) ≈ e atol=1e-10
    end

    @testset "koe2IOE: Equatorial Orbit (i≈0)" begin
        a = 7000.0
        e = 0.1
        i = 1e-12  # Nearly equatorial
        Ω = 1.0
        ω = 0.5
        f = 0.3
        
        koe = SVector{6}(a, e, i, Ω, ω, f)
        IOE = AstroCoords.koe2IOE(koe, μ)
        
        @test IOE[1] ≈ a atol=1e-10
        @test sqrt(IOE[4]^2 + IOE[5]^2) ≈ sin(i/2) atol=1e-14  # Very small
    end

    @testset "koe2IOE: Polar Orbit (i=π/2)" begin
        a = 7000.0
        e = 0.1
        i = π/2
        Ω = 1.0
        ω = 0.5
        f = 0.3
        
        koe = SVector{6}(a, e, i, Ω, ω, f)
        IOE = AstroCoords.koe2IOE(koe, μ)
        
        @test IOE[1] ≈ a atol=1e-10
        @test sqrt(IOE[4]^2 + IOE[5]^2) ≈ sin(π/4) atol=1e-10  # sin(i/2) = sin(π/4)
    end

    @testset "modEqN2IOE: Standard ModEq" begin
        n = sqrt(μ / 7000.0^3)  # Mean motion
        f_meq = 0.01
        g_meq = 0.02
        h_meq = 0.1
        k_meq = 0.05
        L = 1.5
        
        modEq = SVector{6}(n, f_meq, g_meq, h_meq, k_meq, L)
        IOE = AstroCoords.modEqN2IOE(modEq, μ)
        
        @test length(IOE) == 6
        @test IOE[1] ≈ cbrt(μ / n^2) atol=1e-10  # I1 = a
    end

    @testset "IOE2modEqN: Roundtrip Conversion" begin
        # Create IOE state
        I1 = 7000.0
        I2 = 0.01
        I3 = 0.015
        I4 = 0.1
        I5 = 0.05
        I6 = 1.5
        
        IOE_orig = SVector{6}(I1, I2, I3, I4, I5, I6)
        modEq = AstroCoords.IOE2modEqN(IOE_orig, μ)
        IOE_restored = AstroCoords.modEqN2IOE(modEq, μ)
        
        @test IOE_restored[1] ≈ I1 atol=1e-10 rtol=1e-12
        @test IOE_restored[2] ≈ I2 atol=1e-10 rtol=1e-12
        @test IOE_restored[3] ≈ I3 atol=1e-10 rtol=1e-12
        @test IOE_restored[4] ≈ I4 atol=1e-10 rtol=1e-12
        @test IOE_restored[5] ≈ I5 atol=1e-10 rtol=1e-12
        @test IOE_restored[6] ≈ I6 atol=1e-10 rtol=1e-12
    end

    @testset "Type Promotion" begin
        # Test with Int inputs
        koe_int = SVector{6,Int}(7000, 0, 0, 0, 0, 0)
        IOE = AstroCoords.koe2IOE(koe_int, μ)
        @test eltype(IOE) <: AbstractFloat
        
        # Test with Float32 inputs
        koe_f32 = SVector{6,Float32}(7000.0f0, 0.1f0, 0.5f0, 1.0f0, 0.5f0, 0.3f0)
        IOE_f32 = AstroCoords.koe2IOE(koe_f32, Float32(μ))
        @test eltype(IOE_f32) == Float32
    end
end
