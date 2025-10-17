using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

#@lit {citation="Vallado 2013", id="vallado2013", ref="Example 2-5, p.113"}
@testset "cart2koe Extended Coverage" begin
    μ = 398600.4418  # Earth μ [km³/s²]
    
    @testset "Equatorial orbits (i≈0)" begin
        # Test lines 22-48: equatorial_tol parameter
        # Circular equatorial orbit
        r = 7000.0
        v = sqrt(μ / r)
        
        # Exactly equatorial (i = 0)
        state_eq = [r, 0.0, 0.0, 0.0, v, 0.0]
        koe = cart2koe(state_eq, μ; equatorial_tol=1e-15)
        @test koe[1] ≈ r atol=1e-10  # a
        @test koe[2] ≈ 0.0 atol=1e-10  # e
        @test koe[3] ≈ 0.0 atol=1e-15  # i
        @test koe[4] ≈ 0.0 atol=1e-15  # Ω (undefined, should be 0)
        @test koe[5] ≈ 0.0 atol=1e-15  # ω (undefined, should be 0)
        
        # Near-equatorial (within tolerance)
        i_small = 1e-16
        vz = v * sin(i_small)
        vy = v * cos(i_small)
        state_near_eq = [r, 0.0, 0.0, 0.0, vy, vz]
        koe_near = cart2koe(state_near_eq, μ; equatorial_tol=1e-15)
        @test koe_near[3] ≈ 0.0 atol=1e-14  # i treated as 0
        @test koe_near[4] ≈ 0.0 atol=1e-14  # Ω = 0
        
        # Eccentric equatorial orbit
        a = 10000.0
        e = 0.3
        f = π / 4  # true anomaly
        r_mag = a * (1 - e^2) / (1 + e * cos(f))
        v_mag = sqrt(μ * (2/r_mag - 1/a))
        x = r_mag * cos(f)
        y = r_mag * sin(f)
        vx = -v_mag * sin(f) / sqrt(1 - e^2 + 2*e*cos(f) + e^2*cos(f)^2)
        vy = v_mag * (e + cos(f)) / sqrt(1 - e^2 + 2*e*cos(f) + e^2*cos(f)^2)
        state_ecc_eq = [x, y, 0.0, vx, vy, 0.0]
        koe_ecc = cart2koe(state_ecc_eq, μ; equatorial_tol=1e-15)
        @test koe_ecc[2] ≈ e rtol=1e-8  # e
        @test koe_ecc[3] ≈ 0.0 atol=1e-10  # i
    end
    
    @testset "Circular orbits (e≈0)" begin
        # Test lines 51-91: circular_tol parameter
        # Circular inclined orbit
        r = 8000.0
        v = sqrt(μ / r)
        i = π / 6  # 30 degrees
        
        # Position at ascending node
        x = r
        z = 0.0
        vx = 0.0
        vy = v * cos(i)
        vz = v * sin(i)
        state_circ = [x, 0.0, z, vx, vy, vz]
        koe_circ = cart2koe(state_circ, μ; circular_tol=1e-15)
        @test koe_circ[1] ≈ r atol=1e-10  # a
        @test koe_circ[2] ≈ 0.0 atol=1e-10  # e
        @test koe_circ[3] ≈ i rtol=1e-8  # i
        @test koe_circ[5] ≈ 0.0 atol=1e-10  # ω (undefined for circular)
        
        # Near-circular (within tolerance)
        e_small = 1e-16
        a_nc = 8000.0
        r_nc = a_nc * (1 - e_small^2) / (1 + e_small * cos(0))
        v_nc = sqrt(μ * (2/r_nc - 1/a_nc))
        state_near_circ = [r_nc, 0.0, 0.0, 0.0, v_nc, 0.0]
        koe_nc = cart2koe(state_near_circ, μ; circular_tol=1e-15)
        @test koe_nc[2] ≈ 0.0 atol=1e-14  # e treated as 0
    end
    
    @testset "Circular equatorial orbits" begin
        # Test lines 94-133: both tolerances active
        r = 7500.0
        v = sqrt(μ / r)
        state = [r, 0.0, 0.0, 0.0, v, 0.0]
        koe = cart2koe(state, μ; equatorial_tol=1e-15, circular_tol=1e-15)
        @test koe[1] ≈ r atol=1e-10  # a
        @test koe[2] ≈ 0.0 atol=1e-10  # e
        @test koe[3] ≈ 0.0 atol=1e-15  # i
        @test koe[4] ≈ 0.0 atol=1e-15  # Ω
        @test koe[5] ≈ 0.0 atol=1e-15  # ω
        @test 0 <= koe[6] <= 2π  # f
    end
    
    @testset "Near-parabolic orbits (e≈1)" begin
        # Test lines 135-171: numerical stability
        a = 20000.0
        e = 0.99999  # Very close to parabolic
        f = π / 6
        r = a * (1 - e^2) / (1 + e * cos(f))
        v = sqrt(μ * (2/r - 1/a))
        
        # Compute velocity components
        h = sqrt(μ * a * (1 - e^2))
        vr = μ * e * sin(f) / h
        vθ = h / r
        
        x = r * cos(f)
        y = r * sin(f)
        vx = vr * cos(f) - vθ * sin(f)
        vy = vr * sin(f) + vθ * cos(f)
        
        state = [x, y, 0.0, vx, vy, 0.0]
        koe = cart2koe(state, μ)
        @test koe[1] ≈ a rtol=1e-6  # a
        @test koe[2] ≈ e rtol=1e-6  # e
        @test isfinite(koe[6])  # f should be finite
    end
    
    @testset "Hyperbolic orbits (e>1)" begin
        # Test lines 173-210: hyperbolic case
        a = -15000.0  # negative for hyperbolic
        e = 1.5
        f = π / 4
        r = a * (1 - e^2) / (1 + e * cos(f))
        v = sqrt(μ * (2/r - 1/a))
        
        h = sqrt(-μ * a * (1 - e^2))  # Note: a is negative
        vr = μ * e * sin(f) / h
        vθ = h / r
        
        x = r * cos(f)
        y = r * sin(f)
        vx = vr * cos(f) - vθ * sin(f)
        vy = vr * sin(f) + vθ * cos(f)
        
        state = [x, y, 0.0, vx, vy, 0.0]
        koe = cart2koe(state, μ)
        @test koe[1] < 0  # a negative for hyperbolic
        @test koe[2] > 1.0  # e > 1
        @test koe[2] ≈ e rtol=1e-8
    end
    
    @testset "Retrograde orbits (i>90°)" begin
        # Test lines 213-255: retrograde case
        r = 9000.0
        v = sqrt(μ / r)
        i = 2π / 3  # 120 degrees
        
        x = r
        vx = 0.0
        vy = v * cos(i)
        vz = v * sin(i)
        state = [x, 0.0, 0.0, vx, vy, vz]
        koe = cart2koe(state, μ)
        @test koe[3] ≈ i rtol=1e-8  # i
        @test koe[3] > π/2  # retrograde
    end
    
    @testset "Type promotion" begin
        # Test lines 257-287: mixed types
        state_f32 = Float32[7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        μ_f64 = 398600.4418
        koe = cart2koe(state_f32, μ_f64)
        @test eltype(koe) == Float64  # promoted
        
        state_int = [7000, 0, 0, 0, 7, 0]
        koe_int = cart2koe(state_int, μ)
        @test eltype(koe_int) == Float64
    end
    
    @testset "Singularity handlers" begin
        # Test lines 302-339: all edge case branches
        # Already covered in equatorial, circular, etc.
    end
    
    @testset "Round-trip cart→koe→cart" begin
        # Test lines 357-399: verify invertibility
        state_orig = [7000.0, 3000.0, 2000.0, -1.5, 5.0, 4.0]
        koe = cart2koe(state_orig, μ)
        state_back = koe2cart(koe, μ)
        @test state_orig ≈ state_back rtol=1e-10
        
        # Elliptic eccentric
        a = 12000.0
        e = 0.5
        i = π / 4
        Ω = π / 6
        ω = π / 3
        f = π / 2
        koe_test = [a, e, i, Ω, ω, f]
        state = koe2cart(koe_test, μ)
        koe_back = cart2koe(state, μ)
        @test koe_test ≈ koe_back rtol=1e-10
    end
    
    #@lit {citation="Vallado 2013", id="vallado2013", ref="Example 2-5, p.113"}
    @testset "Vallado 2013 Example 2-5" begin
        # Test lines 419-464: published test case
        # Position and velocity from Vallado Example 2-5
        r_vec = [6524.834, 6862.875, 6448.296]  # km
        v_vec = [4.901327, 5.533756, -1.976341]  # km/s
        state = vcat(r_vec, v_vec)
        
        koe = cart2koe(state, μ)
        
        # Expected COE from Vallado (Table 2-4)
        a_exp = 36127.343  # km
        e_exp = 0.832853
        i_exp = deg2rad(87.870)  # radians
        Ω_exp = deg2rad(227.898)  # radians
        ω_exp = deg2rad(53.38)  # radians
        f_exp = deg2rad(92.335)  # radians
        
        # Vallado gives 6 sig figs
        @test koe[1] ≈ a_exp rtol=1e-6
        @test koe[2] ≈ e_exp rtol=1e-6
        @test koe[3] ≈ i_exp rtol=1e-6
        @test koe[4] ≈ Ω_exp rtol=1e-6
        @test koe[5] ≈ ω_exp rtol=1e-6
        @test koe[6] ≈ f_exp rtol=1e-6
    end
    
    @testset "koe→cart edge cases" begin
        # Test lines 472-533: inverse transformation edge cases
        # Circular orbit
        koe_circ = [8000.0, 0.0, π/4, π/6, 0.0, π/3]
        state = koe2cart(koe_circ, μ)
        @test norm(state[1:3]) ≈ 8000.0 rtol=1e-10
        
        # Equatorial orbit
        koe_eq = [10000.0, 0.2, 0.0, 0.0, π/4, π/2]
        state_eq = koe2cart(koe_eq, μ)
        @test abs(state_eq[3]) < 1e-10  # z ≈ 0
        
        # Hyperbolic
        koe_hyp = [-12000.0, 1.3, π/6, π/4, π/3, π/8]
        state_hyp = koe2cart(koe_hyp, μ)
        @test isfinite(norm(state_hyp))
    end
    
    @testset "StaticArrays vs Arrays" begin
        # Test lines 551-595: different array types
        state_vec = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        state_svec = SVector{6}(state_vec)
        
        koe_vec = cart2koe(state_vec, μ)
        koe_svec = cart2koe(state_svec, μ)
        
        @test koe_vec ≈ koe_svec rtol=1e-15
        @test koe_vec isa SVector
        @test koe_svec isa SVector
    end
    
    @testset "Numerical precision limits" begin
        # Test lines 613-692: extreme values
        # Very large orbit
        r_large = 1e8  # km
        v_large = sqrt(μ / r_large)
        state_large = [r_large, 0.0, 0.0, 0.0, v_large, 0.0]
        koe_large = cart2koe(state_large, μ)
        @test isfinite(koe_large[1])
        @test koe_large[1] ≈ r_large rtol=1e-10
        
        # Very small orbit (near collision)
        r_small = 6400.0  # km (Earth radius)
        v_small = sqrt(μ / r_small)
        state_small = [r_small, 0.0, 0.0, 0.0, v_small, 0.0]
        koe_small = cart2koe(state_small, μ)
        @test isfinite(koe_small[1])
        @test koe_small[1] ≈ r_small rtol=1e-10
    end
end
