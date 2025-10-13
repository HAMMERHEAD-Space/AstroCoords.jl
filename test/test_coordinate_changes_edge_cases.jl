using Test
using AstroCoords
using LinearAlgebra

@testset "Coordinate Transformations - Edge Cases" begin
    μ = 398600.4418  # Earth gravitational parameter (km³/s²)

    @testset "cart2koe - Circular Orbit" begin
        # Test with circular orbit (e < circular_tol)
        # Circular orbit at 7000 km altitude
        r = [7000.0, 0.0, 0.0]
        v = [0.0, sqrt(μ/7000.0), 0.0]  # Circular velocity
        
        cart = AstroCoords.Cartesian(r..., v...)
        koe = AstroCoords.cart2koe(cart, μ)
        
        # Check that eccentricity is near zero
        @test koe.e < 1e-10
        
        # Check that semi-major axis matches radius
        @test koe.a ≈ 7000.0 rtol=1e-10
    end

    @testset "cart2koe - Equatorial Orbit" begin
        # Test with equatorial orbit (i < equatorial_tol)
        # Orbit in xy-plane
        a = 8000.0
        e = 0.1
        r_mag = a * (1 - e)  # periapsis
        v_mag = sqrt(μ * (2/r_mag - 1/a))
        
        r = [r_mag, 0.0, 0.0]
        v = [0.0, v_mag, 0.0]
        
        cart = AstroCoords.Cartesian(r..., v...)
        koe = AstroCoords.cart2koe(cart, μ)
        
        # Check that inclination is near zero
        @test koe.i < 1e-10
    end

    @testset "cart2koe - Circular Equatorial Orbit" begin
        # Test with both circular AND equatorial (both singularities)
        r = [7000.0, 0.0, 0.0]
        v = [0.0, sqrt(μ/7000.0), 0.0]
        
        cart = AstroCoords.Cartesian(r..., v...)
        koe = AstroCoords.cart2koe(cart, μ)
        
        # Both e and i should be near zero
        @test koe.e < 1e-10
        @test koe.i < 1e-10
        
        # Other elements should still be defined
        @test isfinite(koe.a)
        @test isfinite(koe.Ω)
    end

    @testset "cart2koe - Retrograde Orbit" begin
        # Test with retrograde orbit (i > 90°)
        a = 8000.0
        e = 0.1
        i = 120.0 * π/180  # 120° inclination
        
        # Build state in perifocal frame then rotate
        r_pf = a * (1 - e)
        v_pf = sqrt(μ * (2/r_pf - 1/a))
        
        # Rotate to include inclination
        r = [r_pf * cos(i), 0.0, r_pf * sin(i)]
        v = [0.0, v_pf * cos(i), v_pf * sin(i)]
        
        cart = AstroCoords.Cartesian(r..., v...)
        koe = AstroCoords.cart2koe(cart, μ)
        
        # Check inclination is > 90°
        @test koe.i > π/2
        @test koe.i ≈ i atol=1e-8
    end

    @testset "koe2cart - Circular Orbit" begin
        # Test converting circular orbit back to Cartesian
        a = 7000.0
        e = 0.0  # exactly circular
        i = 0.5
        Ω = 1.0
        ω = 2.0
        ν = 3.0
        
        koe = AstroCoords.KeplerianOrbitElements(a, e, i, Ω, ω, ν)
        cart = AstroCoords.koe2cart(koe, μ)
        
        # Position magnitude should equal semi-major axis
        r_mag = sqrt(cart.x^2 + cart.y^2 + cart.z^2)
        @test r_mag ≈ a atol=1e-8
        
        # Velocity should be circular velocity
        v_mag = sqrt(cart.ẋ^2 + cart.ẏ^2 + cart.ż^2)
        @test v_mag ≈ sqrt(μ/a) atol=1e-8
    end

    @testset "koe2cart - Equatorial Orbit" begin
        # Test converting equatorial orbit to Cartesian
        a = 7500.0
        e = 0.1
        i = 0.0  # exactly equatorial
        Ω = 1.0
        ω = 2.0
        ν = 3.0
        
        koe = AstroCoords.KeplerianOrbitElements(a, e, i, Ω, ω, ν)
        cart = AstroCoords.koe2cart(koe, μ)
        
        # z-component should be zero (or very small)
        @test abs(cart.z) < 1e-10
        @test abs(cart.ż) < 1e-10
    end

    @testset "koe2cart - Circular Equatorial" begin
        # Test both singularities
        a = 7000.0
        e = 0.0
        i = 0.0
        Ω = 0.0
        ω = 0.0
        ν = π/4  # 45° true anomaly
        
        koe = AstroCoords.KeplerianOrbitElements(a, e, i, Ω, ω, ν)
        cart = AstroCoords.koe2cart(koe, μ)
        
        # Should still produce valid Cartesian state
        @test isfinite(cart.x) && isfinite(cart.y) && isfinite(cart.z)
        @test isfinite(cart.ẋ) && isfinite(cart.ẏ) && isfinite(cart.ż)
    end

    @testset "Round-Trip: cart → koe → cart" begin
        # Test that cart → koe → cart preserves state
        r = [7000.0, 1000.0, 500.0]
        v = [1.0, 7.5, 0.5]
        
        cart_orig = AstroCoords.Cartesian(r..., v...)
        koe = AstroCoords.cart2koe(cart_orig, μ)
        cart_back = AstroCoords.koe2cart(koe, μ)
        
        @test cart_back.x ≈ cart_orig.x rtol=1e-10
        @test cart_back.y ≈ cart_orig.y rtol=1e-10
        @test cart_back.z ≈ cart_orig.z rtol=1e-10
        @test cart_back.ẋ ≈ cart_orig.ẋ rtol=1e-10
        @test cart_back.ẏ ≈ cart_orig.ẏ rtol=1e-10
        @test cart_back.ż ≈ cart_orig.ż rtol=1e-10
    end

    @testset "Energy Conservation" begin
        # Test that energy is preserved in transformations
        a = 8000.0
        e = 0.2
        i = 0.5
        Ω = 1.0
        ω = 2.0
        ν = 3.0
        
        koe = AstroCoords.KeplerianOrbitElements(a, e, i, Ω, ω, ν)
        cart = AstroCoords.koe2cart(koe, μ)
        
        # Compute orbital energy: E = v²/2 - μ/r
        r_mag = sqrt(cart.x^2 + cart.y^2 + cart.z^2)
        v_mag_sq = cart.ẋ^2 + cart.ẏ^2 + cart.ż^2
        energy = v_mag_sq/2 - μ/r_mag
        
        # Energy should match -μ/(2a)
        expected_energy = -μ/(2*a)
        @test energy ≈ expected_energy rtol=1e-10
    end

    @testset "Angular Momentum Conservation" begin
        # Test that angular momentum is preserved
        r = [7000.0, 1000.0, 500.0]
        v = [1.0, 7.5, 0.5]
        
        # Compute angular momentum
        h_orig = cross(r, v)
        h_mag_orig = norm(h_orig)
        
        cart = AstroCoords.Cartesian(r..., v...)
        koe = AstroCoords.cart2koe(cart, μ)
        cart_back = AstroCoords.koe2cart(koe, μ)
        
        r_back = [cart_back.x, cart_back.y, cart_back.z]
        v_back = [cart_back.ẋ, cart_back.ẏ, cart_back.ż]
        h_back = cross(r_back, v_back)
        h_mag_back = norm(h_back)
        
        @test h_mag_back ≈ h_mag_orig rtol=1e-10
    end

    @testset "Anomaly Conversions - Near Circular" begin
        # Test anomaly conversions with e → 0
        e = 1e-10
        M = π/4
        
        E = AstroCoords.meanAnomaly2EccentricAnomaly(M, e)
        ν = AstroCoords.eccentricAnomaly2TrueAnomaly(E, e)
        
        # For circular orbits, M ≈ E ≈ ν
        @test E ≈ M atol=1e-8
        @test ν ≈ M atol=1e-8
    end

    @testset "Anomaly Conversions - Near Parabolic" begin
        # Test with e → 1
        e = 0.9999
        M = 0.1  # Small mean anomaly for parabola
        
        E = AstroCoords.meanAnomaly2EccentricAnomaly(M, e)
        ν = AstroCoords.eccentricAnomaly2TrueAnomaly(E, e)
        
        # Should not have NaN or Inf
        @test isfinite(E)
        @test isfinite(ν)
    end

    @testset "Anomaly Conversions - Hyperbolic" begin
        # Test with e > 1 (hyperbolic)
        e = 1.5
        M = 0.5
        
        # For hyperbolic orbits, use hyperbolic anomaly
        H = AstroCoords.meanAnomaly2HyperbolicAnomaly(M, e)
        ν = AstroCoords.hyperbolicAnomaly2TrueAnomaly(H, e)
        
        @test isfinite(H)
        @test isfinite(ν)
    end
end
