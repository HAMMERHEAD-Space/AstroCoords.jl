using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

#@lit {citation="Solar System Fact Sheet", id="nasa_facts", ref="Earth orbital parameters"}
@testset "Quantity Calculations Coverage" begin
    # Known physical constants
    μ_sun = 1.32712440018e20  # m³/s² (Sun)
    μ_earth = 3.986004418e14  # m³/s² (Earth)
    AU = 1.495978707e11       # m (Astronomical Unit)
    
    @testset "meanMotion - scalar inputs" begin
        # Earth orbit around Sun: a=1AU, period≈365.25 days
        a_earth = AU
        n_earth = meanMotion(a_earth, μ_sun)
        expected_n = 2π / (365.25 * 86400)  # rad/s
        @test n_earth ≈ expected_n atol=1e-10 rtol=1e-6
        
        # LEO: a=6678km, period≈90min
        a_leo = 6.678e6  # m
        n_leo = meanMotion(a_leo, μ_earth)
        expected_period = 90 * 60  # seconds
        @test 2π / n_leo ≈ expected_period atol=10.0 rtol=1e-3
        
        # Type promotion
        n_mixed = meanMotion(Float32(a_leo), μ_earth)
        @test n_mixed isa Float64
    end
    
    @testset "meanMotion - AstroCoord inputs" begin
        # Keplerian elements for Earth orbit
        kep_earth = Keplerian(AU, 0.0167, 0.0, 0.0, 0.0, 0.0)
        n_kep = meanMotion(kep_earth, μ_sun)
        expected_n = 2π / (365.25 * 86400)
        @test n_kep ≈ expected_n atol=1e-10 rtol=1e-6
        
        # Cartesian elements for LEO
        a_leo = 6.678e6
        cart_leo = Cartesian(a_leo, 0.0, 0.0, 0.0, sqrt(μ_earth/a_leo), 0.0)
        n_cart = meanMotion(cart_leo, μ_earth)
        @test n_cart ≈ sqrt(μ_earth / a_leo^3) atol=1e-12 rtol=1e-8
    end
    
    @testset "orbitalPeriod - scalar inputs" begin
        # Earth orbit: period ≈ 365.25 days
        a_earth = AU
        T_earth = orbitalPeriod(a_earth, μ_sun)
        expected_T = 365.25 * 86400  # seconds
        @test T_earth ≈ expected_T atol=100.0 rtol=1e-6
        
        # LEO: period ≈ 90 minutes
        a_leo = 6.678e6
        T_leo = orbitalPeriod(a_leo, μ_earth)
        expected_T_leo = 90 * 60  # seconds
        @test T_leo ≈ expected_T_leo atol=10.0 rtol=1e-3
        
        # Type promotion
        T_mixed = orbitalPeriod(Float32(a_leo), μ_earth)
        @test T_mixed isa Float64
    end
    
    @testset "orbitalPeriod - AstroCoord inputs" begin
        # Keplerian elements
        kep = Keplerian(AU, 0.0167, 0.0, 0.0, 0.0, 0.0)
        T_kep = orbitalPeriod(kep, μ_sun)
        @test T_kep ≈ 365.25 * 86400 atol=100.0 rtol=1e-6
        
        # Modified Equinoctial elements
        p = AU * (1.0 - 0.0167^2)
        modeq = ModEq(p, 0.0, 0.0167, 0.0, 0.0, 0.0)
        T_modeq = orbitalPeriod(modeq, μ_sun)
        @test T_modeq ≈ 365.25 * 86400 atol=100.0 rtol=1e-6
    end
    
    @testset "orbitalNRG - scalar inputs" begin
        # Earth orbit: E = -μ/(2a)
        a_earth = AU
        E_earth = orbitalNRG(a_earth, μ_sun)
        expected_E = -μ_sun / (2.0 * a_earth)
        @test E_earth ≈ expected_E atol=1e-6 rtol=1e-10
        
        # LEO
        a_leo = 6.678e6
        E_leo = orbitalNRG(a_leo, μ_earth)
        expected_E_leo = -μ_earth / (2.0 * a_leo)
        @test E_leo ≈ expected_E_leo atol=1e-6 rtol=1e-10
        
        # Negative energy for bound orbits
        @test E_earth < 0.0
        @test E_leo < 0.0
        
        # Type promotion
        E_mixed = orbitalNRG(Float32(a_leo), μ_earth)
        @test E_mixed isa Float64
    end
    
    @testset "orbitalNRG - AstroCoord inputs" begin
        kep = Keplerian(AU, 0.0167, 0.0, 0.0, 0.0, 0.0)
        E_kep = orbitalNRG(kep, μ_sun)
        @test E_kep ≈ -μ_sun / (2.0 * AU) atol=1e-6 rtol=1e-10
    end
    
    @testset "angularMomentumVector - Cartesian inputs" begin
        # Circular orbit in XY plane
        a = 7.0e6  # m
        v = sqrt(μ_earth / a)
        cart = [a, 0.0, 0.0, 0.0, v, 0.0]
        h_vec = angularMomentumVector(cart)
        
        # Angular momentum should be in Z direction
        @test h_vec[1] ≈ 0.0 atol=1e-6
        @test h_vec[2] ≈ 0.0 atol=1e-6
        @test h_vec[3] > 0.0
        
        # Magnitude check
        expected_h = a * v
        @test norm(h_vec) ≈ expected_h atol=1e-6 rtol=1e-10
        
        # Inclined orbit
        i = π/4  # 45° inclination
        cart_incl = [a, 0.0, 0.0, 0.0, v*cos(i), v*sin(i)]
        h_incl = angularMomentumVector(cart_incl)
        @test norm(h_incl) ≈ expected_h atol=1e-6 rtol=1e-10
    end
    
    @testset "angularMomentumVector - AstroCoord inputs" begin
        # Circular orbit
        a = 7.0e6
        kep = Keplerian(a, 0.0, 0.0, 0.0, 0.0, 0.0)
        h_vec = angularMomentumVector(kep, μ_earth)
        
        expected_h = sqrt(μ_earth * a)
        @test norm(h_vec) ≈ expected_h atol=1e-6 rtol=1e-10
        
        # Elliptic orbit
        kep_ellip = Keplerian(a, 0.2, π/6, 0.0, 0.0, 0.0)
        h_ellip = angularMomentumVector(kep_ellip, μ_earth)
        expected_h_ellip = sqrt(μ_earth * a * (1.0 - 0.2^2))
        @test norm(h_ellip) ≈ expected_h_ellip atol=1e-6 rtol=1e-10
    end
    
    @testset "angularMomentumQuantity - Cartesian inputs" begin
        # Circular orbit
        a = 7.0e6
        v = sqrt(μ_earth / a)
        cart = [a, 0.0, 0.0, 0.0, v, 0.0]
        h_mag = angularMomentumQuantity(cart)
        
        expected_h = a * v
        @test h_mag ≈ expected_h atol=1e-6 rtol=1e-10
        @test h_mag > 0.0
        
        # Should equal norm of vector version
        h_vec = angularMomentumVector(cart)
        @test h_mag ≈ norm(h_vec) atol=1e-12
    end
    
    @testset "angularMomentumQuantity - AstroCoord inputs" begin
        a = 7.0e6
        kep = Keplerian(a, 0.1, π/4, 0.0, 0.0, 0.0)
        h_mag = angularMomentumQuantity(kep, μ_earth)
        
        expected_h = sqrt(μ_earth * a * (1.0 - 0.1^2))
        @test h_mag ≈ expected_h atol=1e-6 rtol=1e-10
        
        # Consistency with vector version
        h_vec = angularMomentumVector(kep, μ_earth)
        @test h_mag ≈ norm(h_vec) atol=1e-12
    end
end
