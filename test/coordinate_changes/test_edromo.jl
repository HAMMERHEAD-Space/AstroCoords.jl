using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "EDromo Transformations" begin
    μ = 398600.4418  # Earth μ [km³/s²]
    
    @testset "cart2EDromo with RegularizedCoordinateConfig and ϕ" begin
        # Test lines 20-24: cart2EDromo with ϕ parameter
        r = 8000.0  # km
        v = sqrt(μ / r)  # circular velocity
        state = [r, 0.0, 0.0, 0.0, v, 0.0]
        
        config = RegularizedCoordinateConfig(state, μ; W=0.0, t₀=0.0, flag_time=LinearTime())
        ϕ = compute_initial_phi(state, μ, config)
        
        edromo_vec = cart2EDromo(state, μ, ϕ, config)
        
        @test length(edromo_vec) == 8
        @test all(isfinite.(edromo_vec))
        @test edromo_vec[3] > 0  # ζ₃ = -1/(2E) > 0 for bound orbits
    end
    
    @testset "DU/TU scaling" begin
        # Test lines 26-36: non-dimensionalization
        state = [7000.0, 3000.0, 2000.0, -1.0, 5.0, 4.0]
        config = RegularizedCoordinateConfig(state, μ; W=0.0)
        ϕ = 0.5
        
        edromo_vec = cart2EDromo(state, μ, ϕ, config)
        
        # Verify state is properly non-dimensionalized
        @test all(isfinite.(edromo_vec))
        
        # Test with different scaling
        config2 = RegularizedCoordinateConfig(
            DU=1e4, TU=100.0, W=0.0, t₀=0.0, flag_time=LinearTime()
        )
        edromo_vec2 = cart2EDromo(state, μ, ϕ, config2)
        @test edromo_vec != edromo_vec2  # Different scaling
    end
    
    @testset "EDromo state vector construction" begin
        # Test lines 42-62: in-plane elements computation
        r = 10000.0
        v = sqrt(μ / r)
        state = [r, 0.0, 0.0, 0.0, v, 0.0]
        config = RegularizedCoordinateConfig(state, μ)
        ϕ = compute_initial_phi(state, μ, config)
        
        edromo_vec = cart2EDromo(state, μ, ϕ, config)
        
        # ζ₁, ζ₂ are in-plane elements
        @test isfinite(edromo_vec[1])  # ζ₁
        @test isfinite(edromo_vec[2])  # ζ₂
        # For circular orbit: ζ₁² + ζ₂² ≈ 0
        @test edromo_vec[1]^2 + edromo_vec[2]^2 < 1e-6
        
        # ζ₃ relates to energy
        @test edromo_vec[3] > 0  # ζ₃ = -1/(2E) > 0
    end
    
    @testset "Time transformation with flag_time" begin
        # Test lines 67-77: different time element formulations
        state = [8000.0, 0.0, 0.0, 0.0, 7.0, 0.0]
        ϕ = 0.3
        
        # PhysicalTime
        config_phys = RegularizedCoordinateConfig(state, μ; t₀=50.0, flag_time=PhysicalTime())
        edromo_phys = cart2EDromo(state, μ, ϕ, config_phys)
        @test edromo_phys[8] ≈ 50.0 / config_phys.TU atol=1e-10  # ζ₈ = t₀/TU
        
        # ConstantTime
        config_const = RegularizedCoordinateConfig(state, μ; t₀=50.0, flag_time=ConstantTime())
        edromo_const = cart2EDromo(state, μ, ϕ, config_const)
        @test isfinite(edromo_const[8])
        
        # LinearTime
        config_lin = RegularizedCoordinateConfig(state, μ; t₀=50.0, flag_time=LinearTime())
        edromo_lin = cart2EDromo(state, μ, ϕ, config_lin)
        @test isfinite(edromo_lin[8])
    end
    
    @testset "Element computation (ζ₁-ζ₈)" begin
        # Test lines 81-117: all eight elements
        state = [9000.0, 4000.0, 3000.0, -2.0, 4.0, 3.0]
        config = RegularizedCoordinateConfig(state, μ)
        ϕ = compute_initial_phi(state, μ, config)
        
        edromo_vec = cart2EDromo(state, μ, ϕ, config)
        
        # Verify all 8 elements are finite
        for i in 1:8
            @test isfinite(edromo_vec[i])
        end
        
        # Quaternion elements (ζ₄-ζ₇) should satisfy normalization
        q_norm = sqrt(edromo_vec[4]^2 + edromo_vec[5]^2 + edromo_vec[6]^2 + edromo_vec[7]^2)
        @test q_norm ≈ 1.0 rtol=1e-10
    end
    
    @testset "Inverse EDromo2cart" begin
        # Test lines 137-189: inverse transformation
        state_orig = [7500.0, 2500.0, 1800.0, -1.5, 5.5, 3.5]
        config = RegularizedCoordinateConfig(state_orig, μ)
        ϕ = compute_initial_phi(state_orig, μ, config)
        
        edromo_vec = cart2EDromo(state_orig, μ, ϕ, config)
        state_back = EDromo2cart(edromo_vec, μ, ϕ, config)
        
        @test state_back[1] ≈ state_orig[1] rtol=1e-10  # x
        @test state_back[2] ≈ state_orig[2] rtol=1e-10  # y
        @test state_back[3] ≈ state_orig[3] rtol=1e-10  # z
        @test state_back[4] ≈ state_orig[4] rtol=1e-10  # vx
        @test state_back[5] ≈ state_orig[5] rtol=1e-10  # vy
        @test state_back[6] ≈ state_orig[6] rtol=1e-10  # vz
    end
    
    @testset "Round-trip verification with tolerance" begin
        # Test complete round-trip transformations
        test_states = [
            [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0],  # circular
            [10000.0, 5000.0, 3000.0, -2.0, 4.0, 3.0],  # elliptic
            [8500.0, 2000.0, 1500.0, -0.8, 6.0, 4.5],  # inclined
        ]
        
        for state in test_states
            config = RegularizedCoordinateConfig(state, μ)
            ϕ = compute_initial_phi(state, μ, config)
            
            edromo_vec = cart2EDromo(state, μ, ϕ, config)
            state_recovered = EDromo2cart(edromo_vec, μ, ϕ, config)
            
            @test state ≈ state_recovered rtol=1e-10
        end
    end
    
    @testset "Numerical stability for near-singular orbits" begin
        # Near-parabolic orbit
        a = 15000.0
        e = 0.99998
        r = a * (1 - e)
        v = sqrt(μ * (2/r - 1/a))
        state_para = [r, 0.0, 0.0, 0.0, v, 0.0]
        config = RegularizedCoordinateConfig(state_para, μ)
        ϕ = compute_initial_phi(state_para, μ, config)
        
        edromo_para = cart2EDromo(state_para, μ, ϕ, config)
        @test all(isfinite.(edromo_para))
        
        # Verify round-trip
        state_back = EDromo2cart(edromo_para, μ, ϕ, config)
        @test state_para ≈ state_back rtol=1e-9
        
        # Circular equatorial (double singularity)
        r_eq = 7500.0
        v_eq = sqrt(μ / r_eq)
        state_eq = [r_eq, 0.0, 0.0, 0.0, v_eq, 0.0]
        config_eq = RegularizedCoordinateConfig(state_eq, μ)
        ϕ_eq = compute_initial_phi(state_eq, μ, config_eq)
        
        edromo_eq = cart2EDromo(state_eq, μ, ϕ_eq, config_eq)
        @test all(isfinite.(edromo_eq))
    end
    
    @testset "Different ϕ (anomaly-like) values" begin
        # Test transformation at different ϕ values
        state = [8000.0, 0.0, 0.0, 0.0, 7.0, 0.0]
        config = RegularizedCoordinateConfig(state, μ)
        
        ϕ_values = [0.0, π/4, π/2, π, 3π/2, 2π]
        
        for ϕ in ϕ_values
            edromo_vec = cart2EDromo(state, μ, ϕ, config)
            @test all(isfinite.(edromo_vec))
            
            # Verify round-trip works at each ϕ
            state_back = EDromo2cart(edromo_vec, μ, ϕ, config)
            @test state ≈ state_back rtol=1e-10
        end
    end
    
    @testset "get_EDromo_time function" begin
        # Test time extraction from EDromo state
        state = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        ϕ = π / 3
        
        # PhysicalTime
        config_phys = RegularizedCoordinateConfig(state, μ; t₀=100.0, flag_time=PhysicalTime())
        edromo_vec = cart2EDromo(state, μ, ϕ, config_phys)
        t_phys = get_EDromo_time(edromo_vec, ϕ, config_phys)
        @test t_phys ≈ 100.0 + edromo_vec[8] * config_phys.TU rtol=1e-10
        
        # ConstantTime
        config_const = RegularizedCoordinateConfig(state, μ; t₀=100.0, flag_time=ConstantTime())
        edromo_vec_const = cart2EDromo(state, μ, ϕ, config_const)
        t_const = get_EDromo_time(edromo_vec_const, ϕ, config_const)
        @test isfinite(t_const)
        
        # LinearTime
        config_lin = RegularizedCoordinateConfig(state, μ; t₀=100.0, flag_time=LinearTime())
        edromo_vec_lin = cart2EDromo(state, μ, ϕ, config_lin)
        t_lin = get_EDromo_time(edromo_vec_lin, ϕ, config_lin)
        @test isfinite(t_lin)
    end
    
    @testset "Type promotion" begin
        # Test with mixed types
        state_f32 = Float32[7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        config = RegularizedCoordinateConfig(state_f32, μ)
        ϕ_f32 = Float32(0.5)
        
        edromo_vec = cart2EDromo(state_f32, μ, ϕ_f32, config)
        @test eltype(edromo_vec) == Float64  # promoted
    end
    
    @testset "W (perturbing potential) effect" begin
        # Test how W affects the transformation
        state = [8000.0, 0.0, 0.0, 0.0, 7.0, 0.0]
        ϕ = 0.5
        
        # Zero W
        config_zero = RegularizedCoordinateConfig(state, μ; W=0.0)
        edromo_zero = cart2EDromo(state, μ, ϕ, config_zero)
        
        # Non-zero W
        W_test = 1e-2
        config_W = RegularizedCoordinateConfig(state, μ; W=W_test)
        edromo_W = cart2EDromo(state, μ, ϕ, config_W)
        
        # Elements should differ when W is non-zero
        @test abs(edromo_zero[3] - edromo_W[3]) > 1e-10  # ζ₃ affected by W
    end
end
