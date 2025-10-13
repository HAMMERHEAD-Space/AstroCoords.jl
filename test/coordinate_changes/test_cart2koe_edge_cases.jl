using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "cart2koe Edge Cases" begin
    μ = 398600.4418  # Earth's gravitational parameter (km^3/s^2)

    @testset "Circular Equatorial Orbit (e≈0, i≈0)" begin
        # Circular orbit in equatorial plane - tests singularity handling at lines 22-48
        r = 7000.0
        v = sqrt(μ / r)
        
        cart = SVector{6}(r, 0.0, 0.0, 0.0, v, 0.0)
        koe = cart2koe(cart, μ)
        
        @test koe[1] ≈ r atol=1e-10  # a = r for circular
        @test koe[2] < 1e-14  # e ≈ 0
        @test koe[3] < 1e-14  # i ≈ 0
        # Ω and ω undefined, f should be well-defined
        @test !isnan(koe[6])
    end

    @testset "Circular Inclined Orbit (e≈0, i>0)" begin
        # Circular orbit with inclination - tests lines 61-91
        r = 7000.0
        v = sqrt(μ / r)
        i = π/4  # 45 degrees
        
        # Position in orbital plane, rotated by inclination
        cart = SVector{6}(r, 0.0, 0.0, 0.0, v*cos(i), v*sin(i))
        koe = cart2koe(cart, μ)
        
        @test koe[1] ≈ r atol=1e-10
        @test koe[2] < 1e-14  # e ≈ 0
        @test koe[3] ≈ i atol=1e-10
        @test !isnan(koe[4])  # Ω should be defined
        # ω undefined for circular orbit
    end

    @testset "Elliptical Equatorial Orbit (e>0, i≈0)" begin
        # Elliptical orbit in equatorial plane - tests lines 112-155
        a = 8000.0
        e = 0.3
        r_peri = a * (1 - e)
        v_peri = sqrt(μ * (2/r_peri - 1/a))
        
        cart = SVector{6}(r_peri, 0.0, 0.0, 0.0, v_peri, 0.0)
        koe = cart2koe(cart, μ)
        
        @test koe[1] ≈ a atol=1e-10
        @test koe[2] ≈ e atol=1e-10
        @test koe[3] < 1e-14  # i ≈ 0
        # Ω undefined, ω and f should be well-defined
        @test !isnan(koe[5])
        @test !isnan(koe[6])
    end

    @testset "Parabolic Trajectory (e=1)" begin
        # Parabolic escape trajectory - tests lines 193-229
        p = 7000.0  # Semi-latus rectum
        e = 1.0
        r = p  # At periapsis
        v = sqrt(2 * μ / r)
        
        cart = SVector{6}(r, 0.0, 0.0, 0.0, v, 0.0)
        koe = cart2koe(cart, μ)
        
        # For parabola, a → ∞, so we check e instead
        @test koe[2] ≈ e atol=1e-10
        @test abs(koe[1]) > 1e6  # Very large a
    end

    @testset "Hyperbolic Trajectory (e>1)" begin
        # Hyperbolic flyby - tests lines 193-229
        a = -7000.0  # Negative for hyperbolic
        e = 1.5
        r_peri = -a * (e - 1)
        v_peri = sqrt(μ * (2/r_peri - 1/a))
        
        cart = SVector{6}(r_peri, 0.0, 0.0, 0.0, v_peri, 0.0)
        koe = cart2koe(cart, μ)
        
        @test koe[1] < 0  # Negative a for hyperbolic
        @test koe[2] > 1.0  # e > 1
        @test koe[2] ≈ e atol=1e-10
    end

    @testset "Retrograde Orbit (i>90°)" begin
        # Retrograde orbit - tests lines 254-287
        r = 7000.0
        v = sqrt(μ / r)
        i = 2*π/3  # 120 degrees (retrograde)
        
        cart = SVector{6}(r, 0.0, 0.0, 0.0, v*cos(i), v*sin(i))
        koe = cart2koe(cart, μ)
        
        @test koe[3] ≈ i atol=1e-10
        @test koe[3] > π/2  # Confirms retrograde
    end

    @testset "Extreme Eccentricity (e=0.99999)" begin
        # Nearly parabolic - tests numerical stability at lines 302-339
        a = 10000.0
        e = 0.99999
        r_peri = a * (1 - e)
        v_peri = sqrt(μ * (2/r_peri - 1/a))
        
        cart = SVector{6}(r_peri, 0.0, 0.0, 0.0, v_peri, 0.0)
        koe = cart2koe(cart, μ)
        
        @test koe[1] ≈ a atol=1e-8  # Relaxed tolerance for extreme case
        @test koe[2] ≈ e atol=1e-8
        @test !isnan(koe[6])  # True anomaly should be computable
    end

    @testset "Custom Tolerances" begin
        # Test with custom equatorial and circular tolerances
        r = 7000.0
        v = sqrt(μ / r)
        
        cart = SVector{6}(r, 0.0, 1e-10, 0.0, v, 0.0)  # Tiny z component
        
        # With default tolerance (1e-15), should detect as equatorial
        koe_default = cart2koe(cart, μ)
        @test koe_default[3] < 1e-14
        
        # With relaxed tolerance, might detect inclination
        koe_relaxed = cart2koe(cart, μ; equatorial_tol=1e-20)
        # Should compute a small but nonzero inclination
        @test koe_relaxed[3] > 0
    end

    @testset "Type Promotion" begin
        # Test with different numeric types
        r = Float32(7000.0)
        v = Float32(sqrt(398600.4418 / 7000.0))
        μ_f32 = Float32(398600.4418)
        
        cart_f32 = SVector{6,Float32}(r, 0.0f0, 0.0f0, 0.0f0, v, 0.0f0)
        koe_f32 = cart2koe(cart_f32, μ_f32)
        
        @test eltype(koe_f32) == Float32
        @test koe_f32[1] ≈ r atol=1e-6  # Float32 precision
    end

    @testset "Zero Angular Momentum Edge Case" begin
        # Radial trajectory (h=0) - extreme edge case
        r = 7000.0
        v_radial = sqrt(2 * μ / r)  # Escape velocity, purely radial
        
        # Position along +x, velocity along +x (radial)
        cart = SVector{6}(r, 0.0, 0.0, v_radial, 0.0, 0.0)
        
        # This should handle gracefully (though orbit is degenerate)
        koe = cart2koe(cart, μ)
        @test !any(isnan, koe)  # No NaNs
    end
end
