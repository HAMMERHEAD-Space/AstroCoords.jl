using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

#@lit {citation="Vallado 2013", id="vallado2013", ref="Example 2-5, p.113"}
@testset "Vallado 2013 — Example 2-5 (p.113): Cartesian to Keplerian" begin
    # Cartesian state from textbook
    r = [6524.834, 6862.875, 6448.296]  # km
    v = [4.901327, 5.533756, -1.976341]  # km/s
    u_cart = vcat(r, v)
    μ = 398600.4418  # km³/s² (Earth)

    # Convert to Keplerian
    u_koe = cart2koe(u_cart, μ)

    # Expected results from textbook
    a_expected = 36127.343  # km
    e_expected = 0.832853
    i_expected = 1.537748  # rad (87.87°)
    Ω_expected = 4.456348  # rad (255.279°)
    ω_expected = 0.634578  # rad (36.345°)
    f_expected = 0.543162  # rad (31.118°)

    # Vallado gives values to 6 significant figures
    @test u_koe[1] ≈ a_expected atol=0.001 rtol=1e-6  # semi-major axis
    @test u_koe[2] ≈ e_expected atol=1e-6 rtol=1e-6  # eccentricity
    @test u_koe[3] ≈ i_expected atol=1e-6 rtol=1e-6  # inclination
    @test u_koe[4] ≈ Ω_expected atol=1e-6 rtol=1e-6  # RAAN
    @test u_koe[5] ≈ ω_expected atol=1e-6 rtol=1e-6  # argument of periapsis
    @test u_koe[6] ≈ f_expected atol=1e-6 rtol=1e-6  # true anomaly
end

#@lit {citation="Vallado 2013", id="vallado2013", ref="Example 2-6, p.118"}
@testset "Vallado 2013 — Example 2-6 (p.118): Keplerian to Cartesian" begin
    # Keplerian elements from textbook
    a = 36127.343  # km
    e = 0.832853
    i = 1.537748  # rad
    Ω = 4.456348  # rad
    ω = 0.634578  # rad
    f = 0.543162  # rad
    u_koe = [a, e, i, Ω, ω, f]
    μ = 398600.4418  # km³/s²

    # Convert to Cartesian
    u_cart = koe2cart(u_koe, μ)

    # Expected results from textbook
    r_expected = [6524.834, 6862.875, 6448.296]  # km
    v_expected = [4.901327, 5.533756, -1.976341]  # km/s

    @test u_cart[1:3] ≈ r_expected atol=0.001 rtol=1e-6
    @test u_cart[4:6] ≈ v_expected atol=0.000001 rtol=1e-6
end

@testset "Edge case: Circular orbit (e=0)" begin
    # Circular orbit test
    μ = 398600.4418
    a = 7000.0  # km
    e = 0.0
    i = 0.5  # rad
    Ω = 1.0  # rad
    ω = 0.0  # rad (undefined for circular, should be set to 0)
    f = 0.5  # rad
    u_koe = [a, e, i, Ω, ω, f]

    # Convert to Cartesian and back
    u_cart = koe2cart(u_koe, μ)
    u_koe_back = cart2koe(u_cart, μ, circular_tol=1e-10)

    # For circular orbits, a and i should be preserved
    @test u_koe_back[1] ≈ a atol=1e-6 rtol=1e-10
    @test u_koe_back[2] ≈ 0.0 atol=1e-12  # eccentricity should be ~0
    @test u_koe_back[3] ≈ i atol=1e-10 rtol=1e-10
    # ω is undefined for circular orbits, typically set to 0
    @test u_koe_back[5] ≈ 0.0 atol=1e-10
end

@testset "Edge case: Equatorial orbit (i=0)" begin
    # Equatorial orbit test
    μ = 398600.4418
    a = 7000.0  # km
    e = 0.1
    i = 0.0
    Ω = 0.0  # rad (undefined for equatorial, should be set to 0)
    ω = 0.5  # rad
    f = 1.0  # rad
    u_koe = [a, e, i, Ω, ω, f]

    # Convert to Cartesian and back
    u_cart = koe2cart(u_koe, μ)
    u_koe_back = cart2koe(u_cart, μ, equatorial_tol=1e-10)

    # For equatorial orbits, a and e should be preserved
    @test u_koe_back[1] ≈ a atol=1e-6 rtol=1e-10
    @test u_koe_back[2] ≈ e atol=1e-10 rtol=1e-10
    @test u_koe_back[3] ≈ 0.0 atol=1e-12  # inclination should be ~0
    # Ω is undefined for equatorial orbits, typically set to 0
    @test u_koe_back[4] ≈ 0.0 atol=1e-10
end

@testset "Edge case: Circular equatorial orbit (e=0, i=0)" begin
    # Circular equatorial orbit
    μ = 398600.4418
    a = 7000.0  # km
    e = 0.0
    i = 0.0
    Ω = 0.0
    ω = 0.0
    f = 0.5  # rad (true longitude in this case)
    u_koe = [a, e, i, Ω, ω, f]

    # Convert to Cartesian and back
    u_cart = koe2cart(u_koe, μ)
    u_koe_back = cart2koe(u_cart, μ, equatorial_tol=1e-10, circular_tol=1e-10)

    # Only a should be well-defined
    @test u_koe_back[1] ≈ a atol=1e-6 rtol=1e-10
    @test u_koe_back[2] ≈ 0.0 atol=1e-12
    @test u_koe_back[3] ≈ 0.0 atol=1e-12
    @test u_koe_back[4] ≈ 0.0 atol=1e-10
    @test u_koe_back[5] ≈ 0.0 atol=1e-10
end

@testset "Edge case: Near-parabolic orbit (e→1)" begin
    # Near-parabolic orbit
    μ = 398600.4418
    a = 20000.0  # km
    e = 0.9999
    i = 0.5  # rad
    Ω = 1.0  # rad
    ω = 0.5  # rad
    f = 0.1  # rad (small true anomaly to stay away from periapsis)
    u_koe = [a, e, i, Ω, ω, f]

    # Convert to Cartesian and back
    u_cart = koe2cart(u_koe, μ)
    u_koe_back = cart2koe(u_cart, μ)

    # Should preserve all elements with reasonable tolerance
    @test u_koe_back[1] ≈ a atol=1e-3 rtol=1e-6
    @test u_koe_back[2] ≈ e atol=1e-8 rtol=1e-6
    @test u_koe_back[3] ≈ i atol=1e-8 rtol=1e-10
    @test u_koe_back[4] ≈ Ω atol=1e-8 rtol=1e-10
    @test u_koe_back[5] ≈ ω atol=1e-8 rtol=1e-10
    @test u_koe_back[6] ≈ f atol=1e-8 rtol=1e-10
end

@testset "Edge case: Retrograde orbit (i>90°)" begin
    # Retrograde orbit
    μ = 398600.4418
    a = 8000.0  # km
    e = 0.2
    i = 2.5  # rad (~143°, retrograde)
    Ω = 1.5  # rad
    ω = 0.8  # rad
    f = 1.2  # rad
    u_koe = [a, e, i, Ω, ω, f]

    # Convert to Cartesian and back
    u_cart = koe2cart(u_koe, μ)
    u_koe_back = cart2koe(u_cart, μ)

    # All elements should be preserved
    @test u_koe_back[1] ≈ a atol=1e-6 rtol=1e-10
    @test u_koe_back[2] ≈ e atol=1e-10 rtol=1e-10
    @test u_koe_back[3] ≈ i atol=1e-10 rtol=1e-10
    @test u_koe_back[4] ≈ Ω atol=1e-10 rtol=1e-10
    @test u_koe_back[5] ≈ ω atol=1e-10 rtol=1e-10
    @test u_koe_back[6] ≈ f atol=1e-10 rtol=1e-10
end

@testset "Hyperbolic orbit (e>1)" begin
    # Hyperbolic trajectory
    μ = 398600.4418
    a = -10000.0  # km (negative for hyperbola)
    e = 1.5
    i = 0.5  # rad
    Ω = 1.0  # rad
    ω = 0.5  # rad
    f = 0.3  # rad (within asymptote limits)
    u_koe = [a, e, i, Ω, ω, f]

    # Convert to Cartesian and back
    u_cart = koe2cart(u_koe, μ)
    u_koe_back = cart2koe(u_cart, μ)

    # All elements should be preserved
    @test u_koe_back[1] ≈ a atol=1e-6 rtol=1e-10
    @test u_koe_back[2] ≈ e atol=1e-10 rtol=1e-10
    @test u_koe_back[3] ≈ i atol=1e-10 rtol=1e-10
    @test u_koe_back[4] ≈ Ω atol=1e-10 rtol=1e-10
    @test u_koe_back[5] ≈ ω atol=1e-10 rtol=1e-10
    @test u_koe_back[6] ≈ f atol=1e-10 rtol=1e-10
end

@testset "Tolerance parameter testing" begin
    μ = 398600.4418
    
    # Test equatorial_tol parameter
    a = 7000.0
    e = 0.1
    i = 1e-16  # Very small inclination
    Ω = 1.0
    ω = 0.5
    f = 1.0
    u_koe = [a, e, i, Ω, ω, f]
    u_cart = koe2cart(u_koe, μ)
    
    # With strict tolerance, should treat as equatorial
    u_koe_strict = cart2koe(u_cart, μ, equatorial_tol=1e-14)
    @test u_koe_strict[3] < 1e-14  # inclination near zero
    @test u_koe_strict[4] ≈ 0.0 atol=1e-10  # RAAN set to zero
    
    # Test circular_tol parameter
    a = 7000.0
    e = 1e-16  # Very small eccentricity
    i = 0.5
    Ω = 1.0
    ω = 0.5
    f = 1.0
    u_koe = [a, e, i, Ω, ω, f]
    u_cart = koe2cart(u_koe, μ)
    
    # With strict tolerance, should treat as circular
    u_koe_strict = cart2koe(u_cart, μ, circular_tol=1e-14)
    @test u_koe_strict[2] < 1e-14  # eccentricity near zero
    @test u_koe_strict[5] ≈ 0.0 atol=1e-10  # AOP set to zero
end

@testset "Coordinate pair conversions: Cartesian ↔ Spherical" begin
    μ = 398600.4418
    
    # Test Cartesian → Spherical → Cartesian
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    u_sphere = cart2sphere(u_cart, μ)
    u_cart_back = sphere2cart(u_sphere, μ)
    
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
end

@testset "Coordinate pair conversions: Cartesian ↔ Cylindrical" begin
    μ = 398600.4418
    
    # Test Cartesian → Cylindrical → Cartesian
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    u_cylind = cart2cylind(u_cart, μ)
    u_cart_back = cylind2cart(u_cylind, μ)
    
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
end

@testset "Coordinate pair conversions: Cartesian ↔ Delaunay" begin
    μ = 398600.4418
    
    # Test Cartesian → Delaunay → Cartesian
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    u_del = cart2delaunay(u_cart, μ)
    u_cart_back = delaunay2cart(u_del, μ)
    
    @test u_cart_back ≈ u_cart atol=1e-6 rtol=1e-8
end

@testset "Coordinate pair conversions: Cartesian ↔ Milankovich" begin
    μ = 398600.4418
    
    # Test Cartesian → Milankovich → Cartesian
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    u_mil = cart2Mil(u_cart, μ)
    u_cart_back = Mil2cart(u_mil, μ)
    
    @test u_cart_back ≈ u_cart atol=1e-6 rtol=1e-8
end

@testset "Coordinate pair conversions: Keplerian ↔ ModEq" begin
    μ = 398600.4418
    
    # Test Keplerian → ModEq → Keplerian
    a = 7000.0
    e = 0.2
    i = 0.5
    Ω = 1.0
    ω = 0.5
    f = 1.0
    u_koe = [a, e, i, Ω, ω, f]
    
    u_modeq = koe2ModEq(u_koe, μ)
    u_koe_back = ModEq2koe(u_modeq, μ)
    
    @test u_koe_back ≈ u_koe atol=1e-8 rtol=1e-10
end

@testset "Coordinate pair conversions: Keplerian ↔ USM7" begin
    μ = 398600.4418
    
    # Test Keplerian → USM7 → Keplerian
    a = 7000.0
    e = 0.2
    i = 0.5
    Ω = 1.0
    ω = 0.5
    f = 1.0
    u_koe = [a, e, i, Ω, ω, f]
    
    u_usm7 = koe2USM7(u_koe, μ)
    u_koe_back = USM72koe(u_usm7, μ)
    
    @test u_koe_back ≈ u_koe atol=1e-8 rtol=1e-10
end

@testset "Coordinate pair conversions: USM7 ↔ USM6" begin
    μ = 398600.4418
    
    # Test USM7 → USM6 → USM7
    # Start from Keplerian
    a = 7000.0
    e = 0.2
    i = 0.5
    Ω = 1.0
    ω = 0.5
    f = 1.0
    u_koe = [a, e, i, Ω, ω, f]
    u_usm7 = koe2USM7(u_koe, μ)
    
    u_usm6 = USM72USM6(u_usm7, μ)
    u_usm7_back = USM62USM7(u_usm6, μ)
    
    @test u_usm7_back ≈ u_usm7 atol=1e-10 rtol=1e-10
end

@testset "Coordinate pair conversions: USM7 ↔ USMEM" begin
    μ = 398600.4418
    
    # Test USM7 → USMEM → USM7
    # Start from Keplerian
    a = 7000.0
    e = 0.2
    i = 0.5
    Ω = 1.0
    ω = 0.5
    f = 1.0
    u_koe = [a, e, i, Ω, ω, f]
    u_usm7 = koe2USM7(u_koe, μ)
    
    u_usmem = USM72USMEM(u_usm7, μ)
    u_usm7_back = USMEM2USM7(u_usmem, μ)
    
    @test u_usm7_back ≈ u_usm7 atol=1e-10 rtol=1e-10
end