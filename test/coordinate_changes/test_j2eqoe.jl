using Test
using AstroCoords
using StaticArrays
using LinearAlgebra

@testset "J2EqOE Transformations" begin
    μ = 398600.4418  # Earth standard gravitational parameter (km^3/s^2)
    
    @testset "koe2IOE - Circular Orbit (e=0)" begin
        # Circular orbit at 7000 km altitude
        a = 7000.0
        e = 0.0
        i = π/4  # 45 degrees
        Ω = π/6  # 30 degrees
        ω = π/3  # 60 degrees
        f = π/4  # 45 degrees
        
        koe = [a, e, i, Ω, ω, f]
        ioe = AstroCoords.koe2IOE(koe, μ)
        
        # For circular orbit, e=0 so I2 and I3 should be zero
        @test ioe[1] ≈ a atol=1e-10
        @test abs(ioe[2]) < 1e-10  # I2 = e*sin(M) ≈ 0
        @test abs(ioe[3]) < 1e-10  # I3 = e*cos(M) ≈ 0
        @test ioe[4] ≈ sin(i/2)*sin(Ω) atol=1e-10
        @test ioe[5] ≈ sin(i/2)*cos(Ω) atol=1e-10
    end
    
    @testset "koe2IOE - Elliptical Orbit" begin
        # Elliptical orbit
        a = 8000.0
        e = 0.3
        i = π/3  # 60 degrees
        Ω = π/4  # 45 degrees
        ω = π/6  # 30 degrees
        f = π/2  # 90 degrees
        
        koe = [a, e, i, Ω, ω, f]
        ioe = AstroCoords.koe2IOE(koe, μ)
        
        @test ioe[1] ≈ a atol=1e-10
        @test norm([ioe[2], ioe[3]]) ≈ e atol=1e-10  # sqrt(I2^2 + I3^2) = e
        @test 2*asin(sqrt(ioe[4]^2 + ioe[5]^2)) ≈ i atol=1e-10
    end
    
    @testset "koe2IOE - Near-Parabolic (e≈1)" begin
        # Near-parabolic orbit
        a = 10000.0
        e = 0.999
        i = π/6
        Ω = π/3
        ω = π/4
        f = π/6
        
        koe = [a, e, i, Ω, ω, f]
        ioe = AstroCoords.koe2IOE(koe, μ)
        
        @test ioe[1] ≈ a atol=1e-10
        @test norm([ioe[2], ioe[3]]) ≈ e atol=1e-8  # Slightly relaxed tolerance for high e
    end
    
    @testset "koe2IOE - Equatorial Orbit (i=0)" begin
        # Equatorial orbit
        a = 7500.0
        e = 0.2
        i = 0.0
        Ω = 0.0  # Undefined for equatorial, use 0
        ω = π/4
        f = π/3
        
        koe = [a, e, i, Ω, ω, f]
        ioe = AstroCoords.koe2IOE(koe, μ)
        
        @test ioe[1] ≈ a atol=1e-10
        @test abs(ioe[4]) < 1e-10  # sin(i/2)*sin(Ω) = 0
        @test abs(ioe[5]) < 1e-10  # sin(i/2)*cos(Ω) = 0
    end
    
    @testset "koe2IOE - Retrograde Orbit (i>90°)" begin
        # Retrograde orbit
        a = 9000.0
        e = 0.15
        i = 2.5  # ~143 degrees
        Ω = π/5
        ω = π/7
        f = π/8
        
        koe = [a, e, i, Ω, ω, f]
        ioe = AstroCoords.koe2IOE(koe, μ)
        
        @test ioe[1] ≈ a atol=1e-10
        @test 2*asin(sqrt(ioe[4]^2 + ioe[5]^2)) ≈ i atol=1e-10
    end
    
    @testset "koe2IOE - Type Promotion" begin
        # Mix Float32 and Float64
        a = Float32(7000.0)
        e = Float64(0.1)
        i = Float32(π/4)
        Ω = Float64(π/6)
        ω = Float32(π/3)
        f = Float64(π/4)
        
        koe = [a, e, i, Ω, ω, f]
        ioe = AstroCoords.koe2IOE(koe, μ)
        
        @test eltype(ioe) == Float64  # Should promote to Float64
        @test ioe[1] ≈ a atol=1e-6
    end
    
    @testset "IOE2koe - Basic Round Trip" begin
        # Test IOE to Keplerian conversion
        a = 8500.0
        e = 0.25
        i = π/3
        Ω = π/4
        ω = π/6
        f = π/2
        
        koe_orig = [a, e, i, Ω, ω, f]
        ioe = AstroCoords.koe2IOE(koe_orig, μ)
        koe_back = AstroCoords.IOE2koe(ioe, μ)
        
        # Round trip should return to original (within tolerance)
        @test koe_back[1] ≈ a atol=1e-10
        @test koe_back[2] ≈ e atol=1e-10
        @test koe_back[3] ≈ i atol=1e-10
        @test koe_back[4] ≈ Ω atol=1e-10
        @test koe_back[5] ≈ ω atol=1e-10
        @test koe_back[6] ≈ f atol=1e-10
    end
    
    @testset "IOE2koe - Equatorial Round Trip" begin
        # Equatorial orbit round trip
        a = 7200.0
        e = 0.1
        i = 0.0
        Ω = 0.0
        ω = π/5
        f = π/4
        
        koe_orig = [a, e, i, Ω, ω, f]
        ioe = AstroCoords.koe2IOE(koe_orig, μ)
        koe_back = AstroCoords.IOE2koe(ioe, μ)
        
        @test koe_back[1] ≈ a atol=1e-10
        @test koe_back[2] ≈ e atol=1e-10
        @test koe_back[3] ≈ i atol=1e-10
    end
    
    @testset "IOE2koe - Circular Round Trip" begin
        # Circular orbit round trip
        a = 6800.0
        e = 0.0
        i = π/6
        Ω = π/4
        ω = π/3
        f = π/2
        
        koe_orig = [a, e, i, Ω, ω, f]
        ioe = AstroCoords.koe2IOE(koe_orig, μ)
        koe_back = AstroCoords.IOE2koe(ioe, μ)
        
        @test koe_back[1] ≈ a atol=1e-10
        @test koe_back[2] ≈ e atol=1e-10
        @test koe_back[3] ≈ i atol=1e-10
    end
    
    @testset "IOE2koe - Near-Parabolic Round Trip" begin
        # Near-parabolic round trip
        a = 12000.0
        e = 0.998
        i = π/8
        Ω = π/7
        ω = π/9
        f = π/12
        
        koe_orig = [a, e, i, Ω, ω, f]
        ioe = AstroCoords.koe2IOE(koe_orig, μ)
        koe_back = AstroCoords.IOE2koe(ioe, μ)
        
        @test koe_back[1] ≈ a atol=1e-8
        @test koe_back[2] ≈ e atol=1e-8
    end
    
    @testset "Edge Case - e=1.0 Boundary" begin
        # Parabolic orbit (e=1.0)
        a = 10000.0
        e = 1.0
        i = π/4
        Ω = π/3
        ω = π/6
        f = π/8
        
        koe = [a, e, i, Ω, ω, f]
        ioe = AstroCoords.koe2IOE(koe, μ)
        
        @test ioe[1] ≈ a atol=1e-10
        @test norm([ioe[2], ioe[3]]) ≈ e atol=1e-10
    end
    
    @testset "Edge Case - i=180° (Polar Retrograde)" begin
        # Polar retrograde orbit
        a = 8000.0
        e = 0.2
        i = π  # 180 degrees
        Ω = π/4
        ω = π/6
        f = π/3
        
        koe = [a, e, i, Ω, ω, f]
        ioe = AstroCoords.koe2IOE(koe, μ)
        koe_back = AstroCoords.IOE2koe(ioe, μ)
        
        @test koe_back[1] ≈ a atol=1e-10
        @test koe_back[2] ≈ e atol=1e-10
        @test koe_back[3] ≈ i atol=1e-10
    end
    
    @testset "Complete Bijection Test" begin
        # Test multiple random cases for bijection
        for _ in 1:10
            a = 6500.0 + 5000.0 * rand()
            e = 0.95 * rand()  # Keep e < 1
            i = π * rand()
            Ω = 2π * rand()
            ω = 2π * rand()
            f = 2π * rand()
            
            koe_orig = [a, e, i, Ω, ω, f]
            ioe = AstroCoords.koe2IOE(koe_orig, μ)
            koe_back = AstroCoords.IOE2koe(ioe, μ)
            
            @test koe_back[1] ≈ koe_orig[1] atol=1e-9
            @test koe_back[2] ≈ koe_orig[2] atol=1e-9
            @test koe_back[3] ≈ koe_orig[3] atol=1e-9
        end
    end
end
