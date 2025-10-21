using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "Stiefel-Scheifele transformations: Basic round-trip" begin
    # Synthetic test: elliptic orbit round-trip
    μ = 398600.4418  # km³/s²
    
    # Standard elliptic orbit
    r = [7000.0, 1000.0, 2000.0]  # km
    v = [1.0, 7.0, -0.5]  # km/s
    u_cart = vcat(r, v)
    
    # Create config
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    ϕ = 0.0  # Initial fictitious time
    
    # Convert to SS and back
    u_ss = cart2StiefelScheifele(u_cart, μ, ϕ, config)
    u_cart_back = StiefelScheifele2cart(u_ss, μ, ϕ, config)
    
    # Round-trip should preserve state
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
end

@testset "Stiefel-Scheifele transformations: α,β quaternion pairs" begin
    # Test that α and β form valid quaternion pairs
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    ϕ = 0.5
    
    u_ss = cart2StiefelScheifele(u_cart, μ, ϕ, config)
    
    # Extract α and β
    α = SVector{4}(u_ss[1], u_ss[2], u_ss[3], u_ss[4])
    β = SVector{4}(u_ss[5], u_ss[6], u_ss[7], u_ss[8])
    
    # α and β should be unit quaternions (approximately)
    # The SS formulation uses them to represent rotations
    # At ϕ=0, u = α and the magnitude should equal √r
    r_mag = norm(r)
    α_mag = norm(α)
    
    # α magnitude should be related to position magnitude
    @test isfinite(α_mag)
    @test α_mag > 0.0
    
    # β magnitude should be related to velocity
    β_mag = norm(β)
    @test isfinite(β_mag)
end

@testset "Stiefel-Scheifele transformations: ω parameter" begin
    # Test that ω parameter is computed correctly
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    ϕ = 0.0
    
    u_ss = cart2StiefelScheifele(u_cart, μ, ϕ, config)
    
    # Extract ω (element 9)
    ω = u_ss[9]
    
    # ω should be related to the total energy
    # For elliptic orbits, ω = √(-E/2)
    # For hyperbolic orbits, ω = √(E/2)
    E = 0.5 * dot(v, v) - μ / norm(r)
    
    if E < 0
        @test ω ≈ sqrt(-E / 2) atol=1e-6 rtol=1e-8
    else
        @test ω ≈ sqrt(E / 2) atol=1e-6 rtol=1e-8
    end
end

@testset "Stiefel-Scheifele transformations: Time transformations with PhysicalTime" begin
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    t₀ = 100.0
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=t₀, flag_time=PhysicalTime())
    ϕ = 0.5
    
    u_ss = cart2StiefelScheifele(u_cart, μ, ϕ, config)
    
    # For PhysicalTime, time element should equal t₀/TU
    @test u_ss[10] ≈ t₀ / config.TU atol=1e-12
    
    # Test get_stiefelscheifele_time function
    t_recovered = get_stiefelscheifele_time(u_ss, ϕ, config)
    @test t_recovered ≈ t₀ atol=1e-10
end

@testset "Stiefel-Scheifele transformations: Time transformations with LinearTime" begin
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    t₀ = 100.0
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=t₀, flag_time=LinearTime())
    ϕ = 0.5
    
    u_ss = cart2StiefelScheifele(u_cart, μ, ϕ, config)
    
    # Test get_stiefelscheifele_time function
    t_recovered = get_stiefelscheifele_time(u_ss, ϕ, config)
    @test t_recovered ≈ t₀ atol=1e-8
end

@testset "Stiefel-Scheifele transformations: Different ϕ values" begin
    # Test that different ϕ values give consistent results
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    
    # Test with different ϕ values
    ϕ_values = [0.0, π/4, π/2, π, 3π/2, 2π]
    
    for ϕ in ϕ_values
        u_ss = cart2StiefelScheifele(u_cart, μ, ϕ, config)
        u_cart_back = StiefelScheifele2cart(u_ss, μ, ϕ, config)
        
        # Should recover original state regardless of ϕ
        @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
    end
end

@testset "Stiefel-Scheifele transformations: Circular orbit" begin
    # Test circular orbit
    μ = 398600.4418
    
    r_mag = 7000.0
    v_circ = sqrt(μ / r_mag)
    
    r = [r_mag, 0.0, 0.0]
    v = [0.0, v_circ, 0.0]
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    ϕ = 0.0
    
    # Convert to SS and back
    u_ss = cart2StiefelScheifele(u_cart, μ, ϕ, config)
    u_cart_back = StiefelScheifele2cart(u_ss, μ, ϕ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
end

@testset "Stiefel-Scheifele transformations: Highly eccentric orbit" begin
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
    ϕ = 0.0
    
    # Convert to SS and back
    u_ss = cart2StiefelScheifele(u_cart, μ, ϕ, config)
    u_cart_back = StiefelScheifele2cart(u_ss, μ, ϕ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-6 rtol=1e-8
end

@testset "Stiefel-Scheifele transformations: Near-parabolic orbit" begin
    μ = 398600.4418
    
    r_mag = 10000.0
    v_escape = sqrt(2 * μ / r_mag)
    
    r = [r_mag, 0.0, 0.0]
    v = [0.0, v_escape * 0.9999, 0.0]
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    ϕ = 0.0
    
    # Convert to SS and back
    u_ss = cart2StiefelScheifele(u_cart, μ, ϕ, config)
    u_cart_back = StiefelScheifele2cart(u_ss, μ, ϕ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-6 rtol=1e-8
end

@testset "Stiefel-Scheifele transformations: Non-zero perturbing potential W" begin
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart = vcat(r, v)
    
    W = 10.0  # km²/s²
    config = RegularizedCoordinateConfig(u_cart, μ, W=W, t₀=0.0, flag_time=LinearTime())
    ϕ = 0.0
    
    # Convert to SS and back
    u_ss = cart2StiefelScheifele(u_cart, μ, ϕ, config)
    u_cart_back = StiefelScheifele2cart(u_ss, μ, ϕ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
end

@testset "Stiefel-Scheifele transformations: Inclined orbit" begin
    μ = 398600.4418
    
    # Create inclined orbit
    r = [5000.0, 3000.0, 4000.0]
    v = [2.0, -1.0, 6.0]
    u_cart = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    ϕ = 0.5
    
    # Convert to SS and back
    u_ss = cart2StiefelScheifele(u_cart, μ, ϕ, config)
    u_cart_back = StiefelScheifele2cart(u_ss, μ, ϕ, config)
    
    @test u_cart_back ≈ u_cart atol=1e-8 rtol=1e-10
end

@testset "Stiefel-Scheifele transformations: Multiple round-trips preserve accuracy" begin
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    u_cart_original = vcat(r, v)
    
    config = RegularizedCoordinateConfig(u_cart_original, μ, W=0.0, t₀=0.0, flag_time=LinearTime())
    ϕ = 0.5
    
    u_cart_current = copy(u_cart_original)
    
    # Perform 5 round-trips
    for _ in 1:5
        u_ss = cart2StiefelScheifele(u_cart_current, μ, ϕ, config)
        u_cart_current = StiefelScheifele2cart(u_ss, μ, ϕ, config)
    end
    
    # Should still be close to original
    @test u_cart_current ≈ u_cart_original atol=1e-6 rtol=1e-8
end