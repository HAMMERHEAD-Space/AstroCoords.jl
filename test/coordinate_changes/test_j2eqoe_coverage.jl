using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "J2 Equinoctial Element Conversions Coverage" begin
    # Test gravitational parameters
    μ_earth = 3.986004418e14  # m³/s²
    
    @testset "koe2IOE - circular orbit (e ≈ 0)" begin
        # Circular orbit at 7000 km altitude
        a = 7.0e6  # m
        e = 0.0
        i = π/6  # 30° inclination
        Ω = π/4  # 45° RAAN
        ω = π/3  # 60° AOP
        f = π/2  # 90° true anomaly
        
        koe = [a, e, i, Ω, ω, f]
        IOE = koe2IOE(koe, μ_earth)
        
        # Check that IOE is returned as SVector
        @test IOE isa SVector{6}
        
        # For circular orbit, e=0 so I2 and I3 should be small
        @test abs(IOE[2]) < 1e-10  # I2 = e*sin(M)
        @test abs(IOE[3]) < 1e-10  # I3 = e*cos(M)
        
        # I1 should equal semi-major axis
        @test IOE[1] ≈ a atol=1e-6
        
        # I4 and I5 involve sin(i/2)
        @test abs(IOE[4]) ≈ sin(i/2) * sin(Ω) atol=1e-10
        @test abs(IOE[5]) ≈ sin(i/2) * cos(Ω) atol=1e-10
    end
    
    @testset "koe2IOE - equatorial orbit (i ≈ 0)" begin
        # Equatorial orbit
        a = 7.0e6
        e = 0.2
        i = 0.0  # Equatorial
        Ω = 0.0
        ω = π/4
        f = π/3
        
        koe = [a, e, i, Ω, ω, f]
        IOE = koe2IOE(koe, μ_earth)
        
        # For equatorial orbit, I4 and I5 should be zero
        @test abs(IOE[4]) < 1e-14  # I4 = sin(i/2)*sin(Ω)
        @test abs(IOE[5]) < 1e-14  # I5 = sin(i/2)*cos(Ω)
        
        # I1 should equal a
        @test IOE[1] ≈ a atol=1e-6
        
        # I2 and I3 should be non-zero (eccentricity present)
        @test abs(IOE[2]) > 0.0
        @test abs(IOE[3]) > 0.0
    end
    
    @testset "koe2IOE - general elliptic orbit" begin
        # General elliptic orbit
        a = 8.0e6
        e = 0.3
        i = π/3  # 60°
        Ω = π/6  # 30°
        ω = π/4  # 45°
        f = π/2  # 90°
        
        koe = [a, e, i, Ω, ω, f]
        IOE = koe2IOE(koe, μ_earth)
        
        # All components should be non-zero
        @test IOE[1] ≈ a atol=1e-6
        @test abs(IOE[2]) > 0.0
        @test abs(IOE[3]) > 0.0
        @test abs(IOE[4]) > 0.0
        @test abs(IOE[5]) > 0.0
        @test IOE[6] > 0.0
    end
    
    @testset "IOE2koe - round-trip conversion" begin
        # Test round-trip: koe → IOE → koe
        a = 7.5e6
        e = 0.15
        i = π/4
        Ω = π/3
        ω = π/6
        f = π/4
        
        koe_orig = [a, e, i, Ω, ω, f]
        IOE = koe2IOE(koe_orig, μ_earth)
        koe_back = IOE2koe(IOE, μ_earth)
        
        # Should recover original Keplerian elements
        @test koe_back[1] ≈ a atol=1e-6 rtol=1e-10  # a
        @test koe_back[2] ≈ e atol=1e-10 rtol=1e-8  # e
        @test koe_back[3] ≈ i atol=1e-10 rtol=1e-8  # i
        @test koe_back[4] ≈ Ω atol=1e-10 rtol=1e-8  # Ω
        @test koe_back[5] ≈ ω atol=1e-10 rtol=1e-8  # ω
        @test koe_back[6] ≈ f atol=1e-10 rtol=1e-8  # f
    end
    
    @testset "IOE2koe - circular orbit recovery" begin
        # Start with circular orbit in IOE
        I1 = 7.0e6  # a
        I2 = 0.0    # e*sin(M)
        I3 = 0.0    # e*cos(M)
        I4 = 0.1    # sin(i/2)*sin(Ω)
        I5 = 0.2    # sin(i/2)*cos(Ω)
        I6 = π/4    # M+ω+Ω
        
        IOE = [I1, I2, I3, I4, I5, I6]
        koe = IOE2koe(IOE, μ_earth)
        
        # Should have zero eccentricity
        @test koe[2] ≈ 0.0 atol=1e-10
        @test koe[1] ≈ I1 atol=1e-6
    end
    
    @testset "IOE2koe - equatorial orbit recovery" begin
        # Equatorial orbit: I4=I5=0
        I1 = 8.0e6
        I2 = 0.05   # Some eccentricity
        I3 = 0.1
        I4 = 0.0    # Equatorial
        I5 = 0.0
        I6 = π/3
        
        IOE = [I1, I2, I3, I4, I5, I6]
        koe = IOE2koe(IOE, μ_earth)
        
        # Should have zero inclination
        @test koe[3] ≈ 0.0 atol=1e-10
        @test koe[1] ≈ I1 atol=1e-6
    end
    
    @testset "modEqN2IOE and IOE2modEqN - round-trip" begin
        # Test ModEq ↔ IOE conversions
        n = sqrt(μ_earth / (7.5e6)^3)
        f = 0.1
        g = 0.2
        h = 0.05
        k = 0.08
        L = π/3
        
        modeqn = [n, f, g, h, k, L]
        IOE = modEqN2IOE(modeqn, μ_earth)
        modeqn_back = IOE2modEqN(IOE, μ_earth)
        
        # Round-trip should preserve values
        @test modeqn_back[1] ≈ n atol=1e-10 rtol=1e-8
        @test modeqn_back[2] ≈ f atol=1e-10 rtol=1e-8
        @test modeqn_back[3] ≈ g atol=1e-10 rtol=1e-8
        @test modeqn_back[4] ≈ h atol=1e-10 rtol=1e-8
        @test modeqn_back[5] ≈ k atol=1e-10 rtol=1e-8
        @test modeqn_back[6] ≈ L atol=1e-10 rtol=1e-8
    end
    
    @testset "Type promotion in conversions" begin
        # Test with mixed types
        koe_mixed = [Float32(7.0e6), Float64(0.1), Float32(π/4), 
                     Float64(π/6), Float32(π/3), Float64(π/2)]
        IOE = koe2IOE(koe_mixed, μ_earth)
        @test IOE isa SVector{6, Float64}
        
        # Float32 μ
        IOE_f32 = koe2IOE(koe_mixed, Float32(μ_earth))
        @test IOE_f32 isa SVector{6, Float64}
    end
    
    @testset "Singularity-free formulation verification" begin
        # Test that IOE elements remain bounded for all orbit types
        
        # Near-circular (e → 0)
        koe_circ = [7.0e6, 1e-10, π/4, π/6, π/3, π/2]
        IOE_circ = koe2IOE(koe_circ, μ_earth)
        @test all(isfinite.(IOE_circ))
        
        # Near-equatorial (i → 0)
        koe_eq = [7.0e6, 0.2, 1e-10, π/6, π/3, π/2]
        IOE_eq = koe2IOE(koe_eq, μ_earth)
        @test all(isfinite.(IOE_eq))
        
        # Both near-singular
        koe_both = [7.0e6, 1e-10, 1e-10, π/6, π/3, π/2]
        IOE_both = koe2IOE(koe_both, μ_earth)
        @test all(isfinite.(IOE_both))
    end
    
    @testset "IOE consistency with mean motion" begin
        # I1 (semi-major axis) should relate to mean motion via Kepler's 3rd law
        a = 7.5e6
        koe = [a, 0.2, π/4, π/6, π/3, π/2]
        IOE = koe2IOE(koe, μ_earth)
        
        n_expected = sqrt(μ_earth / a^3)
        a_from_IOE = IOE[1]
        n_from_IOE = sqrt(μ_earth / a_from_IOE^3)
        
        @test n_from_IOE ≈ n_expected atol=1e-12 rtol=1e-10
    end
    
    @testset "Edge case: zero eccentricity and zero inclination" begin
        # Both e=0 and i=0 (circular equatorial)
        koe = [7.0e6, 0.0, 0.0, 0.0, 0.0, π/4]
        IOE = koe2IOE(koe, μ_earth)
        
        @test IOE[1] ≈ 7.0e6 atol=1e-6
        @test abs(IOE[2]) < 1e-10
        @test abs(IOE[3]) < 1e-10
        @test abs(IOE[4]) < 1e-14
        @test abs(IOE[5]) < 1e-14
        
        # Round-trip
        koe_back = IOE2koe(IOE, μ_earth)
        @test koe_back[2] ≈ 0.0 atol=1e-10  # e
        @test koe_back[3] ≈ 0.0 atol=1e-10  # i
    end
end
