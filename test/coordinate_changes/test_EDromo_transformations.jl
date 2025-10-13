using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "EDromo Transformations" begin
    μ = 398600.4418  # Earth's gravitational parameter (km^3/s^2)

    @testset "cart2EDromo: Standard Elliptical Orbit" begin
        # Standard orbit - tests lines 20-62
        r = 7000.0
        v = sqrt(μ / r)
        cart = SVector{6}(r, 0.0, 0.0, 0.0, v, 0.0)
        
        DU = r
        TU = sqrt(DU^3 / μ)
        W = 0.0
        t₀ = 0.0
        config = RegularizedCoordinateConfig(DU=DU, TU=TU, W=W, t₀=t₀, flag_time=PhysicalTime())
        
        ϕ = compute_initial_phi(cart, μ, config)
        edromo = cart2EDromo(cart, μ, ϕ, config)
        
        @test length(edromo) == 8
        @test !any(isnan, edromo)
        @test edromo[3] > 0  # ζ₃ = -1/(2E) should be positive for elliptical
    end

    @testset "EDromo2cart: Roundtrip Conversion" begin
        # Tests lines 137-189
        r = 8000.0
        v_circ = sqrt(μ / r)
        v = v_circ * 0.95  # Slightly elliptical
        cart_orig = SVector{6}(r, 0.0, 0.0, 0.0, v, v*0.1)
        
        config = RegularizedCoordinateConfig(cart_orig, μ, W=0.0, t₀=0.0, flag_time=PhysicalTime())
        ϕ = compute_initial_phi(cart_orig, μ, config)
        
        edromo = cart2EDromo(cart_orig, μ, ϕ, config)
        cart_restored = EDromo2cart(edromo, μ, ϕ, config)
        
        @test cart_restored[1] ≈ cart_orig[1] atol=1e-10 rtol=1e-12
        @test cart_restored[2] ≈ cart_orig[2] atol=1e-10 rtol=1e-12
        @test cart_restored[3] ≈ cart_orig[3] atol=1e-10 rtol=1e-12
        @test cart_restored[4] ≈ cart_orig[4] atol=1e-10 rtol=1e-12
        @test cart_restored[5] ≈ cart_orig[5] atol=1e-10 rtol=1e-12
        @test cart_restored[6] ≈ cart_orig[6] atol=1e-10 rtol=1e-12
    end

    @testset "Extreme Eccentricities" begin
        # Tests lines 67-85 for numerical edge cases
        
        # Near-circular (e≈0.001)
        r_circ = 7000.0
        v_circ = sqrt(μ / r_circ)
        cart_circ = SVector{6}(r_circ, 0.0, 0.0, 0.0, v_circ*1.001, 0.0)
        
        config_circ = RegularizedCoordinateConfig(cart_circ, μ)
        ϕ_circ = compute_initial_phi(cart_circ, μ, config_circ)
        edromo_circ = cart2EDromo(cart_circ, μ, ϕ_circ, config_circ)
        
        @test !any(isnan, edromo_circ)
        @test sqrt(edromo_circ[1]^2 + edromo_circ[2]^2) < 0.1  # Small eccentricity vector
        
        # High eccentricity (e≈0.999)
        a_ecc = 10000.0
        e_ecc = 0.999
        r_peri = a_ecc * (1 - e_ecc)
        v_peri = sqrt(μ * (2/r_peri - 1/a_ecc))
        cart_ecc = SVector{6}(r_peri, 0.0, 0.0, 0.0, v_peri, 0.0)
        
        config_ecc = RegularizedCoordinateConfig(cart_ecc, μ)
        ϕ_ecc = compute_initial_phi(cart_ecc, μ, config_ecc)
        edromo_ecc = cart2EDromo(cart_ecc, μ, ϕ_ecc, config_ecc)
        
        @test !any(isnan, edromo_ecc)
    end

    @testset "Different ϕ Values" begin
        # Tests lines 101-117 parameter handling
        r = 7500.0
        v = sqrt(μ / r) * 0.9
        cart = SVector{6}(r, 0.0, 0.0, 0.0, v, v*0.2)
        
        config = RegularizedCoordinateConfig(cart, μ)
        
        # Test at different fictitious time values
        for ϕ_test in [0.0, π/4, π/2, π, 3*π/2]
            edromo = cart2EDromo(cart, μ, ϕ_test, config)
            @test !any(isnan, edromo)
            
            # Roundtrip should work at any ϕ
            cart_restored = EDromo2cart(edromo, μ, ϕ_test, config)
            @test cart_restored[1] ≈ cart[1] atol=1e-10 rtol=1e-12
        end
    end

    @testset "Type Promotion and Config Time Handling" begin
        # Tests lines 148-187
        r_f32 = Float32(7000.0)
        v_f32 = Float32(sqrt(398600.4418 / 7000.0))
        cart_f32 = SVector{6,Float32}(r_f32, 0.0f0, 0.0f0, 0.0f0, v_f32, 0.0f0)
        
        DU_f32 = Float32(r_f32)
        TU_f32 = Float32(sqrt(r_f32^3 / 398600.4418))
        config_f32 = RegularizedCoordinateConfig(
            DU=DU_f32, TU=TU_f32, W=0.0f0, t₀=0.0f0, flag_time=PhysicalTime()
        )
        
        ϕ_f32 = Float32(compute_initial_phi(cart_f32, Float32(μ), config_f32))
        edromo_f32 = cart2EDromo(cart_f32, Float32(μ), ϕ_f32, config_f32)
        
        @test eltype(edromo_f32) == Float32
    end

    @testset "Different Time Flags" begin
        # Test different time element formulations
        r = 7000.0
        v = sqrt(μ / r)
        cart = SVector{6}(r, 0.0, 0.0, 0.0, v, 0.0)
        
        # PhysicalTime
        config_phys = RegularizedCoordinateConfig(cart, μ, flag_time=PhysicalTime())
        ϕ = compute_initial_phi(cart, μ, config_phys)
        edromo_phys = cart2EDromo(cart, μ, ϕ, config_phys)
        @test !isnan(edromo_phys[8])  # Time element
        
        # ConstantTime
        config_const = RegularizedCoordinateConfig(cart, μ, flag_time=ConstantTime())
        edromo_const = cart2EDromo(cart, μ, ϕ, config_const)
        @test !isnan(edromo_const[8])
        
        # LinearTime
        config_linear = RegularizedCoordinateConfig(cart, μ, flag_time=LinearTime())
        edromo_linear = cart2EDromo(cart, μ, ϕ, config_linear)
        @test !isnan(edromo_linear[8])
    end

    @testset "get_EDromo_time Function" begin
        # Test time extraction
        r = 7500.0
        v = sqrt(μ / r) * 0.95
        cart = SVector{6}(r, 0.0, 0.0, 0.0, v, v*0.1)
        t₀ = 100.0
        
        config = RegularizedCoordinateConfig(cart, μ, t₀=t₀, flag_time=PhysicalTime())
        ϕ = compute_initial_phi(cart, μ, config)
        edromo = cart2EDromo(cart, μ, ϕ, config)
        
        t = get_EDromo_time(edromo, ϕ, config)
        @test t ≈ t₀ atol=1e-10  # Should return initial time
    end

    @testset "Non-Zero Perturbing Potential" begin
        # Test with W ≠ 0
        r = 7000.0
        v = sqrt(μ / r)
        cart = SVector{6}(r, 0.0, 0.0, 0.0, v, 0.0)
        W = 1000.0  # Some perturbing potential
        
        config = RegularizedCoordinateConfig(cart, μ, W=W)
        ϕ = compute_initial_phi(cart, μ, config)
        edromo = cart2EDromo(cart, μ, ϕ, config)
        
        @test !any(isnan, edromo)
        
        # Roundtrip should work
        cart_restored = EDromo2cart(edromo, μ, ϕ, config)
        @test cart_restored[1] ≈ cart[1] atol=1e-10 rtol=1e-12
    end
end
