using Test
using AstroCoords
using StaticArrays
using LinearAlgebra

@testset "Stiefel-Scheifele Transformations" begin
    # Test parameters
    μ = 398600.4418  # Earth gravitational parameter (km³/s²)
    DU = 6378.137    # km
    TU = 806.8       # seconds
    W = 0.0          # no perturbing potential
    t₀ = 0.0
    
    @testset "RegularizedCoordinateConfig Usage" begin
        # Test config creation with keyword constructor
        config1 = RegularizedCoordinateConfig(
            DU=DU, TU=TU, W=W, t₀=t₀, flag_time=PhysicalTime()
        )
        @test config1.DU ≈ DU
        @test config1.TU ≈ TU
        @test config1.W ≈ W
        @test config1.t₀ ≈ t₀
        @test config1.flag_time isa PhysicalTime
        
        # Test config creation from state
        cart_state = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        config2 = RegularizedCoordinateConfig(cart_state, μ, W=W, t₀=t₀)
        @test config2.DU > 0
        @test config2.TU > 0
        
        # Test different flag_time types
        config_const = RegularizedCoordinateConfig(
            DU=DU, TU=TU, W=W, t₀=t₀, flag_time=ConstantTime()
        )
        @test config_const.flag_time isa ConstantTime
        
        config_linear = RegularizedCoordinateConfig(
            DU=DU, TU=TU, W=W, t₀=t₀, flag_time=LinearTime()
        )
        @test config_linear.flag_time isa LinearTime
    end
    
    @testset "cart2StiefelScheifele - Circular Orbit" begin
        # Circular equatorial orbit
        cart_vec = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        ϕ = 0.0
        config = RegularizedCoordinateConfig(cart_vec, μ)
        
        ss_vec = cart2StiefelScheifele(cart_vec, μ, ϕ, config)
        
        # Check output is SVector{10}
        @test ss_vec isa SVector{10}
        @test length(ss_vec) == 10
        
        # Check α components (first 4 elements)
        @test !isnan(ss_vec[1]) && isfinite(ss_vec[1])
        @test !isnan(ss_vec[2]) && isfinite(ss_vec[2])
        @test !isnan(ss_vec[3]) && isfinite(ss_vec[3])
        @test !isnan(ss_vec[4]) && isfinite(ss_vec[4])
        
        # Check β components (elements 5-8)
        @test !isnan(ss_vec[5]) && isfinite(ss_vec[5])
        @test !isnan(ss_vec[6]) && isfinite(ss_vec[6])
        @test !isnan(ss_vec[7]) && isfinite(ss_vec[7])
        @test !isnan(ss_vec[8]) && isfinite(ss_vec[8])
        
        # Check ω (element 9) - related to energy
        @test ss_vec[9] > 0  # positive for bound orbit
        @test !isnan(ss_vec[9]) && isfinite(ss_vec[9])
        
        # Check time element (element 10)
        @test !isnan(ss_vec[10]) && isfinite(ss_vec[10])
    end
    
    @testset "cart2StiefelScheifele - Elliptical Orbit" begin
        # Elliptical orbit (e=0.3)
        cart_vec = [8000.0, 0.0, 0.0, 0.0, 6.5, 2.0]
        ϕ = 0.0
        config = RegularizedCoordinateConfig(cart_vec, μ)
        
        ss_vec = cart2StiefelScheifele(cart_vec, μ, ϕ, config)
        
        @test ss_vec isa SVector{10}
        @test all(isfinite.(ss_vec))
        @test all(.!isnan.(ss_vec))
        @test ss_vec[9] > 0  # ω positive for elliptical orbit
    end
    
    @testset "cart2StiefelScheifele - Different ϕ Values" begin
        cart_vec = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        config = RegularizedCoordinateConfig(cart_vec, μ)
        
        # Test with ϕ = π/4
        ss1 = cart2StiefelScheifele(cart_vec, μ, π/4, config)
        @test all(isfinite.(ss1))
        
        # Test with ϕ = π/2
        ss2 = cart2StiefelScheifele(cart_vec, μ, π/2, config)
        @test all(isfinite.(ss2))
        
        # Different ϕ should give different states
        @test ss1 != ss2
    end
    
    @testset "StiefelScheifele2cart - Basic Conversion" begin
        # Start with known Cartesian
        cart_vec = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        ϕ = 0.0
        config = RegularizedCoordinateConfig(cart_vec, μ)
        
        # Convert to SS
        ss_vec = cart2StiefelScheifele(cart_vec, μ, ϕ, config)
        
        # Convert back to Cartesian
        cart_back = StiefelScheifele2cart(ss_vec, μ, ϕ, config)
        
        @test cart_back isa SVector{6}
        @test all(isfinite.(cart_back))
    end
    
    @testset "Round-Trip cart→SS→cart" begin
        # Test round-trip conversion accuracy
        cart_orig = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        ϕ = 0.0
        config = RegularizedCoordinateConfig(cart_orig, μ)
        
        ss_vec = cart2StiefelScheifele(cart_orig, μ, ϕ, config)
        cart_back = StiefelScheifele2cart(ss_vec, μ, ϕ, config)
        
        # Check round-trip accuracy with tolerance
        @test cart_back[1] ≈ cart_orig[1] atol=1e-6
        @test cart_back[2] ≈ cart_orig[2] atol=1e-6
        @test cart_back[3] ≈ cart_orig[3] atol=1e-6
        @test cart_back[4] ≈ cart_orig[4] atol=1e-8
        @test cart_back[5] ≈ cart_orig[5] atol=1e-8
        @test cart_back[6] ≈ cart_orig[6] atol=1e-8
    end
    
    @testset "Round-Trip with Elliptical Orbit" begin
        cart_orig = [8000.0, 0.0, 0.0, 0.0, 6.5, 2.0]
        ϕ = π/6
        config = RegularizedCoordinateConfig(cart_orig, μ)
        
        ss_vec = cart2StiefelScheifele(cart_orig, μ, ϕ, config)
        cart_back = StiefelScheifele2cart(ss_vec, μ, ϕ, config)
        
        @test cart_back ≈ cart_orig atol=1e-6 rtol=1e-8
    end
    
    @testset "get_stiefelscheifele_time" begin
        cart_vec = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        ϕ = 0.0
        
        # Test with PhysicalTime
        config_phys = RegularizedCoordinateConfig(
            cart_vec, μ, t₀=100.0, flag_time=PhysicalTime()
        )
        ss_vec_phys = cart2StiefelScheifele(cart_vec, μ, ϕ, config_phys)
        t_phys = get_stiefelscheifele_time(ss_vec_phys, ϕ, config_phys)
        @test t_phys ≈ 100.0 atol=1e-10
        
        # Test with LinearTime
        config_linear = RegularizedCoordinateConfig(
            cart_vec, μ, t₀=100.0, flag_time=LinearTime()
        )
        ss_vec_linear = cart2StiefelScheifele(cart_vec, μ, ϕ, config_linear)
        t_linear = get_stiefelscheifele_time(ss_vec_linear, ϕ, config_linear)
        @test isfinite(t_linear)
    end
    
    @testset "Type Stability" begin
        cart_vec = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        ϕ = 0.0
        config = RegularizedCoordinateConfig(cart_vec, μ)
        
        # Test type promotion
        ss_float64 = cart2StiefelScheifele(cart_vec, μ, ϕ, config)
        @test eltype(ss_float64) === Float64
        
        cart_float32 = Float32.(cart_vec)
        ss_float32 = cart2StiefelScheifele(cart_float32, Float32(μ), Float32(ϕ), config)
        @test eltype(ss_float32) <: AbstractFloat
    end
    
    @testset "Edge Cases" begin
        # Test with perturbing potential
        cart_vec = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        ϕ = 0.0
        config_perturb = RegularizedCoordinateConfig(
            cart_vec, μ, W=10.0, t₀=0.0
        )
        
        ss_perturb = cart2StiefelScheifele(cart_vec, μ, ϕ, config_perturb)
        @test all(isfinite.(ss_perturb))
        
        # Test with non-zero initial time
        config_t0 = RegularizedCoordinateConfig(
            cart_vec, μ, t₀=1000.0
        )
        ss_t0 = cart2StiefelScheifele(cart_vec, μ, ϕ, config_t0)
        @test all(isfinite.(ss_t0))
    end
    
    @testset "Numerical Stability" begin
        # Test near-singularity handling (small radius)
        cart_small = [100.0, 0.0, 0.0, 0.0, 50.0, 0.0]
        ϕ = 0.0
        config_small = RegularizedCoordinateConfig(cart_small, μ)
        
        ss_small = cart2StiefelScheifele(cart_small, μ, ϕ, config_small)
        @test all(isfinite.(ss_small))
        @test all(.!isnan.(ss_small))
    end
end
