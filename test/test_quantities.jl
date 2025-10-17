using Test
using AstroCoords

@testset "Orbital Quantities" begin
    μ = 398600.4418  # Earth's gravitational parameter (km³/s²)

    @testset "Mean Motion" begin
        a = 7000.0  # km
        n = meanMotion(a, μ)
        @test n > 0.0
        @test n ≈ √(μ / a^3) rtol=1e-10

        # Test with Keplerian elements
        kep = Keplerian(a, 0.01, 0.1, 0.2, 0.3, 0.4)
        n_kep = meanMotion(kep, μ)
        @test n_kep ≈ n rtol=1e-10
    end

    @testset "Orbital Period" begin
        a = 7000.0
        T = orbitalPeriod(a, μ)
        @test T > 0.0
        @test T ≈ 2π / √(μ / a^3) rtol=1e-10

        # Test with Keplerian elements
        kep = Keplerian(a, 0.01, 0.1, 0.2, 0.3, 0.4)
        T_kep = orbitalPeriod(kep, μ)
        @test T_kep ≈ T rtol=1e-10
    end

    @testset "orbitalNRG(a, μ) direct computation" begin
        # Test direct computation with positive semi-major axis (elliptic)
        a = 7000.0
        nrg = orbitalNRG(a, μ)
        @test nrg < 0.0  # Elliptic orbits have negative energy
        @test nrg ≈ -μ / (2.0 * a) rtol=1e-10
        
        # Test with larger semi-major axis (less negative energy)
        a_large = 42164.0  # GEO altitude
        nrg_large = orbitalNRG(a_large, μ)
        @test nrg_large < 0.0
        @test nrg_large > nrg  # Less negative (higher energy)
        @test nrg_large ≈ -μ / (2.0 * a_large) rtol=1e-10
    end

    @testset "orbitalNRG from AstroCoord" begin
        # Test with Keplerian elements
        a = 7000.0
        kep = Keplerian(a, 0.01, 0.1, 0.2, 0.3, 0.4)
        nrg_kep = orbitalNRG(kep, μ)
        nrg_direct = orbitalNRG(a, μ)
        @test nrg_kep ≈ nrg_direct rtol=1e-10
        
        # Test with Cartesian coordinates
        # Convert Keplerian to Cartesian
        cart_vec = AstroCoords.koe2cart(AstroCoords.params(kep), μ)
        cart = Cartesian(cart_vec)
        nrg_cart = orbitalNRG(cart, μ)
        @test nrg_cart ≈ nrg_direct rtol=1e-9
    end

    @testset "Energy conservation in transformations" begin
        # Create Keplerian orbit
        a = 7000.0
        e = 0.1
        kep = Keplerian(a, e, 0.5, 1.0, 0.5, 0.3)
        nrg_kep = orbitalNRG(kep, μ)
        
        # Convert to Cartesian and back
        cart = Cartesian(kep, μ)
        nrg_cart = orbitalNRG(cart, μ)
        @test nrg_cart ≈ nrg_kep rtol=1e-9
        
        # Convert back to Keplerian
        kep2 = Keplerian(cart, μ)
        nrg_kep2 = orbitalNRG(kep2, μ)
        @test nrg_kep2 ≈ nrg_kep rtol=1e-9
    end

    @testset "Various semi-major axes (elliptic/hyperbolic)" begin
        # Elliptic orbits (a > 0)
        a_values_elliptic = [6500.0, 7000.0, 8000.0, 10000.0, 42164.0]
        for a in a_values_elliptic
            nrg = orbitalNRG(a, μ)
            @test nrg < 0.0  # Elliptic: negative energy
            @test nrg ≈ -μ / (2.0 * a) rtol=1e-10
        end
        
        # Hyperbolic orbits (a < 0)
        a_values_hyperbolic = [-10000.0, -20000.0, -50000.0]
        for a in a_values_hyperbolic
            nrg = orbitalNRG(a, μ)
            @test nrg > 0.0  # Hyperbolic: positive energy
            @test nrg ≈ -μ / (2.0 * a) rtol=1e-10
        end
    end

    @testset "Type promotion in quantity calculations" begin
        # Test with Float32 and Float64
        a_f32 = 7000.0f0
        μ_f64 = 398600.4418
        
        nrg_mixed = orbitalNRG(a_f32, μ_f64)
        @test nrg_mixed isa Float64  # Should promote to Float64
        
        # Test with all Float32
        μ_f32 = 398600.4418f0
        nrg_f32 = orbitalNRG(a_f32, μ_f32)
        @test nrg_f32 isa Float32
    end

    @testset "Edge cases: a→∞ (escape energy), a→0" begin
        # As a → ∞, energy → 0 from below (escape energy)
        a_large = 1e10
        nrg_large = orbitalNRG(a_large, μ)
        @test nrg_large < 0.0
        @test abs(nrg_large) < 1e-4  # Very close to zero
        
        # As a → 0, energy → -∞
        a_small = 1e-3
        nrg_small = orbitalNRG(a_small, μ)
        @test nrg_small < -1e8  # Very negative
        
        # Verify relationship: smaller a means more negative energy
        a1 = 7000.0
        a2 = 14000.0
        nrg1 = orbitalNRG(a1, μ)
        nrg2 = orbitalNRG(a2, μ)
        @test nrg1 < nrg2  # Smaller orbit has more negative energy
    end

    @testset "Angular Momentum Vector" begin
        # Test with Cartesian state
        r = [7000.0, 0.0, 0.0]
        v = [0.0, 7.5, 0.0]
        u = vcat(r, v)
        h_vec = angularMomentumVector(u)
        @test length(h_vec) == 3
        @test h_vec[3] ≈ 7000.0 * 7.5 rtol=1e-10  # h = r × v, z-component

        # Test with AstroCoord
        cart = Cartesian(u)
        h_vec_coord = angularMomentumVector(cart, μ)
        @test h_vec_coord ≈ h_vec rtol=1e-10
    end

    @testset "Angular Momentum Quantity" begin
        # Test with circular orbit in xy-plane
        r = [7000.0, 0.0, 0.0]
        v = [0.0, √(μ/7000.0), 0.0]  # Circular velocity
        u = vcat(r, v)
        h = angularMomentumQuantity(u)
        @test h > 0.0
        @test h ≈ 7000.0 * √(μ/7000.0) rtol=1e-9

        # Test with AstroCoord
        cart = Cartesian(u)
        h_coord = angularMomentumQuantity(cart, μ)
        @test h_coord ≈ h rtol=1e-9
    end

    @testset "Vis-viva equation verification" begin
        # Verify vis-viva: v² = μ(2/r - 1/a)
        a = 7000.0
        e = 0.1
        f = π/4  # True anomaly
        
        kep = Keplerian(a, e, 0.5, 1.0, 0.5, f)
        cart = Cartesian(kep, μ)
        
        # Extract position and velocity
        r_vec = [cart.x, cart.y, cart.z]
        v_vec = [cart.ẋ, cart.ẏ, cart.ż]
        
        r_mag = √(cart.x^2 + cart.y^2 + cart.z^2)
        v_mag = √(cart.ẋ^2 + cart.ẏ^2 + cart.ż^2)
        
        # Vis-viva equation
        v_squared_visviva = μ * (2.0/r_mag - 1.0/a)
        v_squared_actual = v_mag^2
        
        @test v_squared_actual ≈ v_squared_visviva rtol=1e-9
        
        # Verify energy from vis-viva
        nrg_visviva = v_squared_actual/2.0 - μ/r_mag
        nrg_direct = orbitalNRG(a, μ)
        @test nrg_visviva ≈ nrg_direct rtol=1e-9
    end
end
