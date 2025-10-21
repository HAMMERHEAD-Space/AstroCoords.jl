using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "EDromo transformations: Basic round-trip (elliptic)" begin
    # Synthetic test: elliptic orbit round-trip
    μ = 398600.4418  # km³/s²
    
    # Standard elliptic orbit
    r = [7000.0, 1000.0, 2000.0]  # km
    v = [1.0, 7.0, -0.5]  # km/s
    u_cart = vcat(r, v)
    
    # Create config and compute initial phi
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    ϕ = compute_initial_phi(u_cart, μ, config)
    
    # Convert to EDromo and back
    u_edromo = cart2EDromo(u_cart, μ, ϕ, config)
    u_cart_back = EDromo2cart(u_edromo, μ, ϕ, config)
    
    # Round-trip should preserve state
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
end

@testset "EDromo transformations: Circular orbit" begin
    # Test circular orbit
    μ = 398600.4418
    
    # Circular orbit (e ≈ 0)
    r_mag = 7000.0
    v_circ = sqrt(μ / r_mag)
    
    r = [r_mag, 0.0, 0.0]
    v = [0.0, v_circ, 0.0]
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    ϕ = compute_initial_phi(u_cart, μ, config)
    
    # Convert to EDromo and back
    u_edromo = cart2EDromo(u_cart, μ, ϕ, config)
    u_cart_back = EDromo2cart(u_edromo, μ, ϕ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
    
    # For circular orbit, ζ₁² + ζ₂² should be small
    @test sqrt(u_edromo[1]^2 + u_edromo[2]^2) < 0.1
end

@testset "EDromo transformations: Near-parabolic orbit" begin
    # Test near-parabolic orbit (E ≈ 0)
    μ = 398600.4418
    
    r_mag = 10000.0
    v_escape = sqrt(2 * μ / r_mag)
    
    r = [r_mag, 0.0, 0.0]
    v = [0.0, v_escape * 0.9999, 0.0]  # Very close to escape velocity
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    ϕ = compute_initial_phi(u_cart, μ, config)
    
    # Convert to EDromo and back
    u_edromo = cart2EDromo(u_cart, μ, ϕ, config)
    u_cart_back = EDromo2cart(u_edromo, μ, ϕ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-6 rtol=1e-8
end

@testset "EDromo transformations: ϕ parameter effects" begin
    # Test that different ϕ values give consistent results
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    
    # Test with different ϕ values
    ϕ_values = [0.0, π/4, π/2, π, 3π/2]
    
    for ϕ in ϕ_values
        u_edromo = cart2EDromo(u_cart, μ, ϕ, config)
        u_cart_back = EDromo2cart(u_edromo, μ, ϕ, config)
        
        # Should recover original state regardless of ϕ
        @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
    end
end

@testset "EDromo transformations: Time element with PhysicalTime" begin
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    t₀ = 100.0
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=t₀, flag_time=PhysicalTime())
    ϕ = compute_initial_phi(u_cart, μ, config)
    
    u_edromo = cart2EDromo(u_cart, μ, ϕ, config)
    
    # For PhysicalTime, ζ₈ should equal t₀/TU
    @test u_edromo[8] ≈ t₀ / config.TU atol=1e-12
    
    # Test get_EDromo_time function
    t_recovered = get_EDromo_time(u_edromo, ϕ, config)
    @test t_recovered ≈ t₀ atol=1e-10
end

@testset "EDromo transformations: Time element with ConstantTime" begin
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    t₀ = 100.0
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=t₀, flag_time=ConstantTime())
    ϕ = compute_initial_phi(u_cart, μ, config)
    
    u_edromo = cart2EDromo(u_cart, μ, ϕ, config)
    
    # Test get_EDromo_time function
    t_recovered = get_EDromo_time(u_edromo, ϕ, config)
    @test t_recovered ≈ t₀ atol=1e-8
end

@testset "EDromo transformations: Time element with LinearTime" begin
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    t₀ = 100.0
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=t₀, flag_time=LinearTime())
    ϕ = compute_initial_phi(u_cart, μ, config)
    
    u_edromo = cart2EDromo(u_cart, μ, ϕ, config)
    
    # Test get_EDromo_time function
    t_recovered = get_EDromo_time(u_edromo, ϕ, config)
    @test t_recovered ≈ t₀ atol=1e-8
end

@testset "EDromo transformations: Non-zero perturbing potential W" begin
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    W = 10.0  # km²/s²
    config = RegularizedCoordinateConfig(u_cart, μ, W=W, t₀=0.0, flag_time=LinearTime())
    ϕ = compute_initial_phi(u_cart, μ, config)
    
    # Convert to EDromo and back
    u_edromo = cart2EDromo(u_cart, μ, ϕ, config)
    u_cart_back = EDromo2cart(u_edromo, μ, ϕ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
end

@testset "EDromo transformations: RegularizedCoordinateConfig with various DU/TU" begin
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    # Test with different DU/TU scales
    DU_scales = [1000.0, 5000.0, 10000.0]
    
    for DU in DU_scales
        TU = sqrt(DU^3 / μ)
        config = RegularizedCoordinateConfig(DU=DU, TU=TU, W=0.0, t₀=0.0, flag_time=LinearTime())
        ϕ = compute_initial_phi(u_cart, μ, config)
        
        # Convert to EDromo and back
        u_edromo = cart2EDromo(u_cart, μ, ϕ, config)
        u_cart_back = EDromo2cart(u_edromo, μ, ϕ, config)
        
        # Result should be independent of scaling
        @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
    end
end

@testset "EDromo transformations: Highly eccentric orbit" begin
    μ = 398600.4418
    
    # Highly eccentric orbit near periapsis
    r_p = 7000.0
    e = 0.95
    a = r_p / (1 - e)
    v_p = sqrt(μ * (1 + e) / r_p)
    
    r = [r_p, 0.0, 0.0]
    v = [0.0, v_p, 0.0]
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    ϕ = compute_initial_phi(u_cart, μ, config)
    
    # Convert to EDromo and back
    u_edromo = cart2EDromo(u_cart, μ, ϕ, config)
    u_cart_back = EDromo2cart(u_edromo, μ, ϕ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-6 rtol=1e-8
end

@testset "EDromo transformations: Inclined orbit" begin
    # Test orbit with significant inclination
    μ = 398600.4418
    
    # Create inclined orbit
    r = [5000.0, 3000.0, 4000.0]
    v = [2.0, -1.0, 6.0]
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    ϕ = compute_initial_phi(u_cart, μ, config)
    
    # Convert to EDromo and back
    u_edromo = cart2EDromo(u_cart, μ, ϕ, config)
    u_cart_back = EDromo2cart(u_edromo, μ, ϕ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
end

@testset "EDromo transformations: compute_initial_phi consistency" begin
    # Test that compute_initial_phi produces reasonable values
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    ϕ = compute_initial_phi(u_cart, μ, config)
    
    # ϕ should be finite and real
    @test isfinite(ϕ)
    @test !isnan(ϕ)
    
    # ϕ should be in a reasonable range (typically -2π to 2π)
    @test abs(ϕ) < 10.0
end

@testset "EDromo transformations: Multiple round-trips preserve accuracy" begin
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart_original = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart_original, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    
    u_cart_current = copy(u_cart_original)
    
    # Perform 5 round-trips
    for _ in 1:5
        ϕ = compute_initial_phi(u_cart_current, μ, config)
        u_edromo = cart2EDromo(u_cart_current, μ, ϕ, config)
        u_cart_current = EDromo2cart(u_edromo, μ, ϕ, config)
    end
    
    # Should still be close to original
    @test u_cart_current ≈ u_cart_original atol=1e-6 rtol=1e-8
end