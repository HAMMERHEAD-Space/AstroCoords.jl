using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "Kustaanheimo-Stiefel Transformations" begin
    μ = 398600.4418  # Earth μ [km³/s²]
    
    @testset "cart2KS with RegularizedCoordinateConfig" begin
        # Test lines 16-20: cart2KS with config
        r = 7000.0  # km
        v = sqrt(μ / r)  # circular velocity
        state = [r, 0.0, 0.0, 0.0, v, 0.0]
        
        config = RegularizedCoordinateConfig(state, μ; W=0.0, t₀=0.0, flag_time=LinearTime())
        ks_vec = cart2KS(state, μ, config)
        
        @test length(ks_vec) == 10
        @test all(isfinite.(ks_vec))
        # Verify KS position magnitude relates to physical r
        u_mag_sq = ks_vec[1]^2 + ks_vec[2]^2 + ks_vec[3]^2 + ks_vec[4]^2
        @test u_mag_sq * config.DU ≈ r rtol=1e-10
    end
    
    @testset "W computation and scaling" begin
        # Test lines 25-33: W parameter in transformation
        r = 8000.0
        v = sqrt(μ / r)
        state = [r, 0.0, 0.0, 0.0, v, 0.0]
        
        # Test with zero W
        config_zero = RegularizedCoordinateConfig(state, μ; W=0.0)
        ks_zero = cart2KS(state, μ, config_zero)
        @test isfinite(ks_zero[9])  # h (energy)
        
        # Test with non-zero W
        W_test = 1e-3
        config_W = RegularizedCoordinateConfig(state, μ; W=W_test)
        ks_W = cart2KS(state, μ, config_W)
        @test isfinite(ks_W[9])
        # Energy should differ when W is non-zero
        @test abs(ks_zero[9] - ks_W[9]) > 1e-10
    end
    
    @testset "KS matrix transformations" begin
        # Test lines 39-48: KS position computation
        state = [6000.0, 3000.0, 2000.0, -1.0, 5.0, 3.0]
        config = RegularizedCoordinateConfig(state, μ)
        ks_vec = cart2KS(state, μ, config)
        
        # Verify KS coordinates follow bilinear relations
        u1, u2, u3, u4 = ks_vec[1:4]
        r_x = (u1^2 - u2^2 - u3^2 + u4^2) * config.DU
        r_y = 2 * (u1*u2 - u3*u4) * config.DU
        r_z = 2 * (u1*u3 + u2*u4) * config.DU
        
        @test r_x ≈ state[1] rtol=1e-10
        @test r_y ≈ state[2] rtol=1e-10
        @test r_z ≈ state[3] rtol=1e-10
    end
    
    @testset "Velocity transformations" begin
        # Test lines 52-64: KS velocity computation
        state = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        config = RegularizedCoordinateConfig(state, μ)
        ks_vec = cart2KS(state, μ, config)
        
        # KS velocity components should be finite
        @test all(isfinite.(ks_vec[5:8]))
        
        # Test velocity magnitude relationship
        u = SVector{4}(ks_vec[1:4])
        u_dot = SVector{4}(ks_vec[5:8])
        @test isfinite(dot(u, u_dot))
    end
    
    @testset "Time-element computation" begin
        # Test lines 67-87: time element with different flag_time
        state = [8000.0, 0.0, 0.0, 0.0, 7.0, 0.0]
        
        # PhysicalTime
        config_phys = RegularizedCoordinateConfig(state, μ; flag_time=PhysicalTime())
        ks_phys = cart2KS(state, μ, config_phys)
        @test ks_phys[10] ≈ 0.0 atol=1e-15  # τ = t₀/TU
        
        # LinearTime
        config_lin = RegularizedCoordinateConfig(state, μ; flag_time=LinearTime())
        ks_lin = cart2KS(state, μ, config_lin)
        @test isfinite(ks_lin[10])
    end
    
    @testset "Inverse transformation KS2cart" begin
        # Test lines 83-99: KS2cart function
        state_orig = [7500.0, 2000.0, 1500.0, -0.5, 6.0, 4.0]
        config = RegularizedCoordinateConfig(state_orig, μ)
        
        ks_vec = cart2KS(state_orig, μ, config)
        state_back = KS2cart(ks_vec, μ, config)
        
        @test state_back[1] ≈ state_orig[1] rtol=1e-10  # x
        @test state_back[2] ≈ state_orig[2] rtol=1e-10  # y
        @test state_back[3] ≈ state_orig[3] rtol=1e-10  # z
        @test state_back[4] ≈ state_orig[4] rtol=1e-10  # vx
        @test state_back[5] ≈ state_orig[5] rtol=1e-10  # vy
        @test state_back[6] ≈ state_orig[6] rtol=1e-10  # vz
    end
    
    @testset "Round-trip cart→KS→cart verification" begin
        # Test lines 102-106: complete round-trip
        states = [
            [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0],  # circular
            [10000.0, 5000.0, 3000.0, -2.0, 4.0, 3.0],  # elliptic
            [8000.0, 0.0, 0.0, 0.0, sqrt(μ/8000.0), 0.0],  # circular
        ]
        
        for state in states
            config = RegularizedCoordinateConfig(state, μ)
            ks_vec = cart2KS(state, μ, config)
            state_recovered = KS2cart(ks_vec, μ, config)
            @test state ≈ state_recovered rtol=1e-10
        end
    end
    
    @testset "Numerical stability near collision (r→0)" begin
        # Test numerical behavior at small r
        r_small = 6400.0  # Near Earth surface
        v_small = sqrt(μ / r_small)
        state_small = [r_small, 0.0, 0.0, 0.0, v_small, 0.0]
        
        config = RegularizedCoordinateConfig(state_small, μ)
        ks_vec = cart2KS(state_small, μ, config)
        @test all(isfinite.(ks_vec))
        
        # Verify round-trip still works
        state_back = KS2cart(ks_vec, μ, config)
        @test state_small ≈ state_back rtol=1e-9
    end
    
    @testset "Different orbit types" begin
        # Elliptic
        a = 12000.0
        e = 0.5
        r = a * (1 - e)
        v = sqrt(μ * (2/r - 1/a))
        state_elliptic = [r, 0.0, 0.0, 0.0, v, 0.0]
        config_e = RegularizedCoordinateConfig(state_elliptic, μ)
        ks_e = cart2KS(state_elliptic, μ, config_e)
        @test ks_e[9] < 0  # h < 0 for bound orbits
        
        # Hyperbolic
        a_hyp = -15000.0
        e_hyp = 1.5
        r_hyp = abs(a_hyp) * (e_hyp - 1)
        v_hyp = sqrt(μ * (2/r_hyp - 1/a_hyp))
        state_hyperbolic = [r_hyp, 0.0, 0.0, 0.0, v_hyp, 0.0]
        config_h = RegularizedCoordinateConfig(state_hyperbolic, μ)
        ks_h = cart2KS(state_hyperbolic, μ, config_h)
        @test ks_h[9] > 0  # h > 0 for unbound orbits
        
        # Near-parabolic
        a_para = 20000.0
        e_para = 0.99999
        r_para = a_para * (1 - e_para)
        v_para = sqrt(μ * (2/r_para - 1/a_para))
        state_parabolic = [r_para, 0.0, 0.0, 0.0, v_para, 0.0]
        config_p = RegularizedCoordinateConfig(state_parabolic, μ)
        ks_p = cart2KS(state_parabolic, μ, config_p)
        @test isfinite(ks_p[9])
    end
    
    @testset "get_KS_time function" begin
        # Test time extraction from KS state
        state = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        
        # PhysicalTime
        config_phys = RegularizedCoordinateConfig(state, μ; t₀=100.0, flag_time=PhysicalTime())
        ks_vec = cart2KS(state, μ, config_phys)
        t_phys = get_KS_time(ks_vec, config_phys)
        @test t_phys ≈ 100.0 rtol=1e-10
        
        # LinearTime
        config_lin = RegularizedCoordinateConfig(state, μ; t₀=100.0, flag_time=LinearTime())
        ks_vec_lin = cart2KS(state, μ, config_lin)
        t_lin = get_KS_time(ks_vec_lin, config_lin)
        @test isfinite(t_lin)
        @test abs(t_lin - 100.0) < 10.0  # Should be near t₀
    end
    
    @testset "Type promotion" begin
        # Test with mixed types
        state_f32 = Float32[7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        config = RegularizedCoordinateConfig(state_f32, μ)
        ks_vec = cart2KS(state_f32, μ, config)
        @test eltype(ks_vec) == Float64  # promoted
    end
end
