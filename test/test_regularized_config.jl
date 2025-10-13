using Test
using AstroCoords
using Dates

@testset "RegularizedCoordinateConfig" begin
    @testset "Full Constructor" begin
        # Test full constructor with all parameters
        TT = 100.0
        DU = 1000.0
        TT2 = DateTime(2024, 1, 1)
        W = 1.5
        
        config = AstroCoords.RegularizedCoordinateConfig(TT, DU, TT2, W)
        
        @test config.TT == TT
        @test config.DU == DU
        @test config.TT2 == TT2
        @test config.W == W
    end

    @testset "Convenience Constructor from State" begin
        # Test constructor that computes DU/TU from state
        μ = 398600.4418  # km³/s² (Earth)
        
        # Create test state (circular orbit at ~7000 km)
        r = [7000.0, 0.0, 0.0]  # position (km)
        v = [0.0, 7.5, 0.0]     # velocity (km/s)
        
        TT2 = 0.0
        W = 1.0
        
        config = AstroCoords.RegularizedCoordinateConfig(r, v, μ, TT2, W)
        
        # DU should be magnitude of r
        @test config.DU ≈ 7000.0 atol=1e-10
        
        # TU should be computed from orbital period
        # For circular orbit: period = 2π√(a³/μ), so TU ~ period/2π
        expected_TU = sqrt(7000.0^3 / μ)
        @test config.TT ≈ expected_TU atol=1e-8
        
        @test config.TT2 == TT2
        @test config.W == W
    end

    @testset "Different Time Types" begin
        # Test with DateTime
        config1 = AstroCoords.RegularizedCoordinateConfig(
            100.0, 1000.0, DateTime(2024, 1, 1), 1.0
        )
        @test config1.TT2 isa DateTime
        
        # Test with Float
        config2 = AstroCoords.RegularizedCoordinateConfig(
            100.0, 1000.0, 0.0, 1.0
        )
        @test config2.TT2 isa Float64
        
        # Test with Integer
        config3 = AstroCoords.RegularizedCoordinateConfig(
            100.0, 1000.0, 0, 1.0
        )
        @test config3.TT2 isa Int
    end

    @testset "Type Promotion" begin
        # Test type promotion in constructor
        config = AstroCoords.RegularizedCoordinateConfig(
            Float32(100.0), Float64(1000.0), 0.0, Float32(1.0)
        )
        
        # Should promote to common type
        @test config.TT isa Float64 || config.TT isa Float32
        @test config.DU isa Float64
        @test config.W isa Float64 || config.W isa Float32
    end

    @testset "Field Access" begin
        # Test that all fields are accessible
        config = AstroCoords.RegularizedCoordinateConfig(
            100.0, 1000.0, 0.0, 1.5
        )
        
        @test config.TT == 100.0
        @test config.DU == 1000.0
        @test config.TT2 == 0.0
        @test config.W == 1.5
    end

    @testset "Config with Elliptical Orbit" begin
        # Test with elliptical orbit parameters
        μ = 398600.4418
        
        # Elliptical orbit (a=8000, e=0.1)
        r = [7200.0, 0.0, 0.0]
        v = [0.0, 7.8, 0.0]
        
        config = AstroCoords.RegularizedCoordinateConfig(r, v, μ, 0.0, 1.0)
        
        @test config.DU > 0
        @test config.TT > 0
        @test isfinite(config.DU)
        @test isfinite(config.TT)
    end

    @testset "Config with High Eccentricity" begin
        # Test with highly eccentric orbit
        μ = 398600.4418
        
        # Eccentric orbit at apogee (low velocity)
        r = [20000.0, 0.0, 0.0]
        v = [0.0, 3.0, 0.0]
        
        config = AstroCoords.RegularizedCoordinateConfig(r, v, μ, 0.0, 1.0)
        
        @test config.DU ≈ 20000.0 atol=1e-10
        @test config.TT > 0
        @test isfinite(config.TT)
    end

    @testset "Config Consistency" begin
        # Test that manually created and state-derived configs are consistent
        μ = 398600.4418
        r = [7000.0, 0.0, 0.0]
        v = [0.0, 7.5, 0.0]
        
        # Create from state
        config1 = AstroCoords.RegularizedCoordinateConfig(r, v, μ, 0.0, 1.0)
        
        # Create manually with computed values
        DU = norm(r)
        TU = sqrt(DU^3 / μ)
        config2 = AstroCoords.RegularizedCoordinateConfig(TU, DU, 0.0, 1.0)
        
        @test config1.TT ≈ config2.TT atol=1e-10
        @test config1.DU ≈ config2.DU atol=1e-10
    end
end
