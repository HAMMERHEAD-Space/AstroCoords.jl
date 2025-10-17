using Test
using AstroCoords
using StaticArrays

@testset "RegularizedCoordinateConfig" begin
    @testset "Constructor with explicit parameters" begin
        # Test lines 19-22: explicit parameter constructor
        config = RegularizedCoordinateConfig(
            DU=1e7, TU=1e3, W=0.0, t₀=0.0, flag_time=PhysicalTime()
        )
        @test config.DU ≈ 1e7 atol=1e-10
        @test config.TU ≈ 1e3 atol=1e-10
        @test config.W ≈ 0.0 atol=1e-15
        @test config.t₀ ≈ 0.0 atol=1e-15
        @test config.flag_time isa PhysicalTime
    end

    @testset "Constructor from state vector" begin
        # Test lines 38-46: computing DU/TU from state
        μ = 398600.4418  # Earth μ [km³/s²]
        
        # Circular orbit at 7000 km
        r = 7000.0  # km
        v = sqrt(μ / r)  # circular velocity
        state = [r, 0.0, 0.0, 0.0, v, 0.0]
        
        config = RegularizedCoordinateConfig(state, μ; W=0.0, t₀=0.0, flag_time=LinearTime())
        
        # Verify DU is computed from position magnitude
        @test config.DU ≈ r atol=1e-10
        # Verify TU follows Kepler's third law: TU = sqrt(DU³/μ)
        expected_TU = sqrt(r^3 / μ)
        @test config.TU ≈ expected_TU rtol=1e-10
        @test config.W ≈ 0.0 atol=1e-15
        @test config.t₀ ≈ 0.0 atol=1e-15
        @test config.flag_time isa LinearTime
    end

    @testset "Different AbstractTimeType implementations" begin
        # Test lines 65-68: different flag_time types
        state = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        μ = 398600.4418
        
        config_physical = RegularizedCoordinateConfig(state, μ; flag_time=PhysicalTime())
        @test config_physical.flag_time isa PhysicalTime
        
        config_constant = RegularizedCoordinateConfig(state, μ; flag_time=ConstantTime())
        @test config_constant.flag_time isa ConstantTime
        
        config_linear = RegularizedCoordinateConfig(state, μ; flag_time=LinearTime())
        @test config_linear.flag_time isa LinearTime
    end

    @testset "Type promotion" begin
        # Test lines 88-92: mixing number types
        config_mixed = RegularizedCoordinateConfig(
            DU=1e7,      # Float64
            TU=1000,     # Int64
            W=0.0f0,     # Float32
            t₀=0,        # Int64
            flag_time=PhysicalTime()
        )
        
        # All fields should be promoted to Float64
        @test config_mixed.DU isa Float64
        @test config_mixed.TU isa Float64
        @test config_mixed.W isa Float64
        @test config_mixed.t₀ isa Float64
    end

    @testset "W (perturbing potential) computation" begin
        # Test lines 94-97: W parameter accuracy
        μ = 398600.4418
        state = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        
        # Test with non-zero W (perturbing potential)
        W_test = 1e-5  # Small perturbation
        config = RegularizedCoordinateConfig(state, μ; W=W_test)
        @test config.W ≈ W_test atol=1e-15
        
        # Test with zero W (no perturbation)
        config_zero = RegularizedCoordinateConfig(state, μ; W=0.0)
        @test config_zero.W ≈ 0.0 atol=1e-15
        
        # Test with negative W (valid for certain perturbations)
        W_neg = -1e-6
        config_neg = RegularizedCoordinateConfig(state, μ; W=W_neg)
        @test config_neg.W ≈ W_neg atol=1e-15
    end

    @testset "flag_time boolean behavior" begin
        # Test lines 99-100: flag_time usage
        state = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        μ = 398600.4418
        
        # PhysicalTime: τ = t directly
        config_phys = RegularizedCoordinateConfig(state, μ; flag_time=PhysicalTime())
        @test config_phys.flag_time isa PhysicalTime
        
        # LinearTime: τ includes additional terms
        config_lin = RegularizedCoordinateConfig(state, μ; flag_time=LinearTime())
        @test config_lin.flag_time isa LinearTime
        
        # ConstantTime: τ computed with constant offset
        config_const = RegularizedCoordinateConfig(state, μ; flag_time=ConstantTime())
        @test config_const.flag_time isa ConstantTime
    end

    @testset "Edge cases" begin
        μ = 398600.4418
        
        # Test lines 102-104: near-parabolic orbit (e ≈ 1)
        a = 10000.0
        e = 0.9999  # Near-parabolic
        r = a * (1 - e)  # periapsis
        v = sqrt(μ * (2/r - 1/a))  # vis-viva
        state_parabolic = [r, 0.0, 0.0, 0.0, v, 0.0]
        config_para = RegularizedCoordinateConfig(state_parabolic, μ)
        @test config_para.DU ≈ r atol=1e-10
        @test isfinite(config_para.TU)
        
        # Hyperbolic orbit (e > 1)
        a_hyp = -10000.0  # negative for hyperbolic
        e_hyp = 1.5
        r_hyp = abs(a_hyp) * (e_hyp - 1)  # periapsis
        v_hyp = sqrt(μ * (2/r_hyp - 1/a_hyp))
        state_hyperbolic = [r_hyp, 0.0, 0.0, 0.0, v_hyp, 0.0]
        config_hyp = RegularizedCoordinateConfig(state_hyperbolic, μ)
        @test config_hyp.DU ≈ r_hyp atol=1e-10
        @test isfinite(config_hyp.TU)
        
        # Very small orbit (LEO)
        r_leo = 6571.0  # km (200 km altitude)
        v_leo = sqrt(μ / r_leo)
        state_leo = [r_leo, 0.0, 0.0, 0.0, v_leo, 0.0]
        config_leo = RegularizedCoordinateConfig(state_leo, μ)
        @test config_leo.DU ≈ r_leo atol=1e-10
        @test config_leo.TU ≈ sqrt(r_leo^3 / μ) rtol=1e-10
    end

    @testset "compute_characteristic_scales" begin
        # Test the helper function directly
        μ = 398600.4418
        r_mag = 8000.0
        state = [r_mag, 0.0, 0.0, 0.0, 7.0, 0.0]
        
        DU, TU = compute_characteristic_scales(state, μ)
        @test DU ≈ r_mag atol=1e-10
        @test TU ≈ sqrt(r_mag^3 / μ) rtol=1e-10
        
        # Test with non-axis-aligned position
        state_3d = [5000.0, 3000.0, 4000.0, 1.0, 2.0, 3.0]
        r_3d = sqrt(5000^2 + 3000^2 + 4000^2)
        DU_3d, TU_3d = compute_characteristic_scales(state_3d, μ)
        @test DU_3d ≈ r_3d rtol=1e-10
        @test TU_3d ≈ sqrt(r_3d^3 / μ) rtol=1e-10
    end
end
