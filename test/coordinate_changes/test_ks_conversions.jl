using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "KS transformations: Basic round-trip" begin
    # Synthetic test: round-trip conversion
    μ = 398600.4418  # km³/s²
    
    # Standard elliptic orbit
    r = [7000.0, 1000.0, 2000.0]  # km
    v = [1.0, 7.0, -0.5]  # km/s
    u_cart = vcat(r, v)
    
    # Create config
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    
    # Convert to KS and back
    u_ks = cart2KS(u_cart, μ, config)
    u_cart_back = KS2cart(u_ks, μ, config)
    
    # Round-trip should preserve state within numerical precision
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
end

@testset "KS transformations: Periapsis passage" begin
    # Test at periapsis (r·v = 0)
    μ = 398600.4418
    
    # Position at periapsis (velocity perpendicular to position)
    r = [6578.137, 0.0, 0.0]  # km (near Earth surface)
    v = [0.0, 7.784, 0.0]  # km/s (perpendicular)
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    
    # Convert to KS and back
    u_ks = cart2KS(u_cart, μ, config)
    u_cart_back = KS2cart(u_ks, μ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
    
    # At periapsis, r·v should be near zero
    @test dot(u_cart_back[1:3], u_cart_back[4:6]) ≈ 0.0 atol=1e-6
end

@testset "KS transformations: Zero energy orbit (parabolic)" begin
    # Test near-zero energy orbit
    μ = 398600.4418
    
    # Construct near-parabolic orbit (v² ≈ 2μ/r)
    r_mag = 10000.0  # km
    v_escape = sqrt(2 * μ / r_mag)
    
    r = [r_mag, 0.0, 0.0]
    v = [0.0, v_escape * 0.9999, 0.0]  # Very close to escape velocity
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    
    # Convert to KS and back
    u_ks = cart2KS(u_cart, μ, config)
    u_cart_back = KS2cart(u_ks, μ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-6 rtol=1e-8
    
    # Energy should be near zero
    E = 0.5 * dot(v, v) - μ / r_mag
    @test abs(E) < 100.0  # Near zero energy (small compared to μ/r)
end

@testset "KS transformations: Time element transformations" begin
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    # Test PhysicalTime
    t₀ = 100.0
    config_phys = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=t₀, flag_time=PhysicalTime())
    u_ks_phys = cart2KS(u_cart, μ, config_phys)
    
    # For PhysicalTime, τ should equal t₀/TU
    @test u_ks_phys[10] ≈ t₀ / config_phys.TU atol=1e-12
    
    # Test get_KS_time function
    t_recovered = get_KS_time(u_ks_phys, config_phys)
    @test t_recovered ≈ t₀ atol=1e-10
    
    # Test LinearTime
    config_linear = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=t₀, flag_time=LinearTime())
    u_ks_linear = cart2KS(u_cart, μ, config_linear)
    
    # For LinearTime, τ includes additional term
    t_recovered_linear = get_KS_time(u_ks_linear, config_linear)
    @test t_recovered_linear ≈ t₀ atol=1e-8
end

@testset "KS transformations: RegularizedCoordinateConfig with various DU/TU" begin
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    # Test with different DU/TU scales
    DU_scales = [1000.0, 5000.0, 10000.0]
    
    for DU in DU_scales
        TU = sqrt(DU^3 / μ)
        config = RegularizedCoordinateConfig(DU=DU, TU=TU, W=0.0, t₀=0.0, flag_time=LinearTime())
        
        # Convert to KS and back
        u_ks = cart2KS(u_cart, μ, config)
        u_cart_back = KS2cart(u_ks, μ, config)
        
        # Result should be independent of scaling
        @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
    end
end

@testset "KS transformations: Non-zero perturbing potential W" begin
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    # Test with non-zero perturbing potential
    W = 10.0  # km²/s²
    config = RegularizedCoordinateConfig(u_cart, μ, W=W, t₀=0.0, flag_time=LinearTime())
    
    # Convert to KS and back
    u_ks = cart2KS(u_cart, μ, config)
    u_cart_back = KS2cart(u_ks, μ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
    
    # The energy element h should account for W
    E_total = 0.5 * dot(v, v) - μ / norm(r) + W
    @test u_ks[9] ≈ -E_total atol=1e-6 rtol=1e-8
end

@testset "KS transformations: Highly eccentric orbit" begin
    # Test numerical stability with high eccentricity
    μ = 398600.4418
    
    # Highly eccentric orbit near periapsis
    r_p = 7000.0  # km (periapsis radius)
    e = 0.95
    a = r_p / (1 - e)
    v_p = sqrt(μ * (1 + e) / r_p)  # Velocity at periapsis
    
    r = [r_p, 0.0, 0.0]
    v = [0.0, v_p, 0.0]
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    
    # Convert to KS and back
    u_ks = cart2KS(u_cart, μ, config)
    u_cart_back = KS2cart(u_ks, μ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-6 rtol=1e-8
end

@testset "KS transformations: Edge case x[1] < 0" begin
    # Test the branch where x[1] < 0 in cart2KS
    μ = 398600.4418
    
    # Position with negative x-component
    r = [-7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    
    # Convert to KS and back
    u_ks = cart2KS(u_cart, μ, config)
    u_cart_back = KS2cart(u_ks, μ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
end

@testset "KS transformations: Hyperbolic trajectory" begin
    # Test hyperbolic orbit (E > 0)
    μ = 398600.4418
    
    # Hyperbolic trajectory
    r_mag = 10000.0
    v_mag = sqrt(3 * μ / r_mag)  # v > v_escape
    
    r = [r_mag, 0.0, 0.0]
    v = [0.0, v_mag, 0.0]
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    
    # Convert to KS and back
    u_ks = cart2KS(u_cart, μ, config)
    u_cart_back = KS2cart(u_ks, μ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
    
    # Energy should be positive
    E = 0.5 * v_mag^2 - μ / r_mag
    @test E > 0.0
    @test u_ks[9] ≈ -E atol=1e-6 rtol=1e-8
end

@testset "KS transformations: Multiple round-trips preserve accuracy" begin
    # Test that multiple conversions don't accumulate significant errors
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart_original = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart_original, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    
    u_cart_current = copy(u_cart_original)
    
    # Perform 5 round-trips
    for _ in 1:5
        u_ks = cart2KS(u_cart_current, μ, config)
        u_cart_current = KS2cart(u_ks, μ, config)
    end
    
    # Should still be close to original
    @test u_cart_current ≈ u_cart_original atol=1e-6 rtol=1e-8
end