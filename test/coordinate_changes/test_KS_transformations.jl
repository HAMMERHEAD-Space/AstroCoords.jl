using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "Kustaanheimo-Stiefel Transformations" begin
    μ = 398600.4418  # Earth's gravitational parameter (km^3/s^2)

    @testset "cart2KS: Standard Orbit with Config" begin
        # Tests lines 16-33 - setup and matrix operations
        r = 7000.0
        v = sqrt(μ / r)
        cart = SVector{6}(r, 0.0, 0.0, 0.0, v, 0.0)
        
        config = RegularizedCoordinateConfig(cart, μ, W=0.0, t₀=0.0, flag_time=PhysicalTime())
        ks = cart2KS(cart, μ, config)
        
        @test length(ks) == 10
        @test !any(isnan, ks)
        
        # Check that KS position gives correct radius
        r_ks = ks[1]^2 + ks[2]^2 + ks[3]^2 + ks[4]^2
        @test r_ks ≈ (r / config.DU) atol=1e-12
    end

    @testset "KS2cart: Roundtrip Validation" begin
        # Tests lines 99-106 - roundtrip with strict tolerances
        r = 8000.0
        v_circ = sqrt(μ / r)
        v = v_circ * 0.9
        cart_orig = SVector{6}(r, 0.0, 0.0, 0.0, v, v*0.15)
        
        config = RegularizedCoordinateConfig(cart_orig, μ)
        ks = cart2KS(cart_orig, μ, config)
        cart_restored = KS2cart(ks, μ, config)
        
        @test cart_restored[1] ≈ cart_orig[1] atol=1e-10 rtol=1e-12
        @test cart_restored[2] ≈ cart_orig[2] atol=1e-10 rtol=1e-12
        @test cart_restored[3] ≈ cart_orig[3] atol=1e-10 rtol=1e-12
        @test cart_restored[4] ≈ cart_orig[4] atol=1e-10 rtol=1e-12
        @test cart_restored[5] ≈ cart_orig[5] atol=1e-10 rtol=1e-12
        @test cart_restored[6] ≈ cart_orig[6] atol=1e-10 rtol=1e-12
    end

    @testset "Edge Cases: Near-Circular Orbit" begin
        # Tests lines 39-67 for L-matrix stability with r≈a
        r = 7000.0
        v = sqrt(μ / r) * 1.001  # Slightly more than circular velocity
        cart = SVector{6}(r, 0.0, 0.0, 0.0, v, 0.0)
        
        config = RegularizedCoordinateConfig(cart, μ)
        ks = cart2KS(cart, μ, config)
        
        @test !any(isnan, ks)
        
        cart_restored = KS2cart(ks, μ, config)
        @test cart_restored[1] ≈ cart[1] atol=1e-10 rtol=1e-12
    end

    @testset "Edge Cases: Highly Elliptical at Periapsis" begin
        # Tests lines 39-67 - r << a
        a = 10000.0
        e = 0.9
        r_peri = a * (1 - e)
        v_peri = sqrt(μ * (2/r_peri - 1/a))
        cart = SVector{6}(r_peri, 0.0, 0.0, 0.0, v_peri, 0.0)
        
        config = RegularizedCoordinateConfig(cart, μ)
        ks = cart2KS(cart, μ, config)
        
        @test !any(isnan, ks)
        
        cart_restored = KS2cart(ks, μ, config)
        @test cart_restored[1] ≈ cart[1] atol=1e-9  # Slightly relaxed for extreme case
    end

    @testset "Edge Cases: Highly Elliptical at Apoapsis" begin
        # Tests lines 39-67 - r >> a at apoapsis
        a = 10000.0
        e = 0.9
        r_apo = a * (1 + e)
        v_apo = sqrt(μ * (2/r_apo - 1/a))
        cart = SVector{6}(r_apo, 0.0, 0.0, 0.0, v_apo, 0.0)
        
        config = RegularizedCoordinateConfig(cart, μ)
        ks = cart2KS(cart, μ, config)
        
        @test !any(isnan, ks)
        
        cart_restored = KS2cart(ks, μ, config)
        @test cart_restored[1] ≈ cart[1] atol=1e-9
    end

    @testset "Different W Values" begin
        # Tests lines 83-96 - regularization parameter
        r = 7500.0
        v = sqrt(μ / r) * 0.95
        cart = SVector{6}(r, 0.0, 0.0, 0.0, v, v*0.1)
        
        for W_test in [0.0, 500.0, 1000.0, 2000.0]
            config = RegularizedCoordinateConfig(cart, μ, W=W_test)
            ks = cart2KS(cart, μ, config)
            
            @test !any(isnan, ks)
            
            # Roundtrip should work with any W
            cart_restored = KS2cart(ks, μ, config)
            @test cart_restored[1] ≈ cart[1] atol=1e-10 rtol=1e-12
        end
    end

    @testset "Time Conversion Flags" begin
        # Ensures proper s↔t handling
        r = 7000.0
        v = sqrt(μ / r)
        cart = SVector{6}(r, 0.0, 0.0, 0.0, v, 0.0)
        t₀ = 100.0
        
        # PhysicalTime
        config_phys = RegularizedCoordinateConfig(
            cart, μ, t₀=t₀, flag_time=PhysicalTime()
        )
        ks_phys = cart2KS(cart, μ, config_phys)
        t_phys = get_KS_time(ks_phys, config_phys)
        @test t_phys ≈ t₀ atol=1e-10
        
        # LinearTime
        config_linear = RegularizedCoordinateConfig(
            cart, μ, t₀=t₀, flag_time=LinearTime()
        )
        ks_linear = cart2KS(cart, μ, config_linear)
        t_linear = get_KS_time(ks_linear, config_linear)
        @test !isnan(t_linear)
    end

    @testset "get_KS_time Function" begin
        r = 7500.0
        v = sqrt(μ / r) * 0.95
        cart = SVector{6}(r, 0.0, 0.0, 0.0, v, v*0.1)
        t₀ = 50.0
        
        config = RegularizedCoordinateConfig(cart, μ, t₀=t₀, flag_time=PhysicalTime())
        ks = cart2KS(cart, μ, config)
        
        t = get_KS_time(ks, config)
        @test t ≈ t₀ atol=1e-10
    end

    @testset "Negative x-Position Branch" begin
        # Tests the x < 0 branch in cart2KS (lines 44-49)
        r = 7000.0
        v = sqrt(μ / r)
        cart_neg = SVector{6}(-r, 0.0, 0.0, 0.0, v, 0.0)
        
        config = RegularizedCoordinateConfig(cart_neg, μ)
        ks = cart2KS(cart_neg, μ, config)
        
        @test !any(isnan, ks)
        
        cart_restored = KS2cart(ks, μ, config)
        @test cart_restored[1] ≈ -r atol=1e-10 rtol=1e-12
    end

    @testset "Type Promotion" begin
        r_f32 = Float32(7000.0)
        v_f32 = Float32(sqrt(398600.4418 / 7000.0))
        cart_f32 = SVector{6,Float32}(r_f32, 0.0f0, 0.0f0, 0.0f0, v_f32, 0.0f0)
        μ_f32 = Float32(μ)
        
        config_f32 = RegularizedCoordinateConfig(
            cart_f32, μ_f32, W=0.0f0, t₀=0.0f0, flag_time=PhysicalTime()
        )
        ks_f32 = cart2KS(cart_f32, μ_f32, config_f32)
        
        @test eltype(ks_f32) == Float32
    end
end
