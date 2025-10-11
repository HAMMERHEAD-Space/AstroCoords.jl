using Test
using AstroCoords
using LinearAlgebra

#@lit {citation="Vallado 2013", id="vallado2013", ref="Example 2-7, Algorithm 9"}
@testset "J2 Equinoctial Orbital Elements" begin
    μ = 398600.4418  # Earth gravitational parameter (km³/s²)
    
    @testset "Classical to J2EqOE conversions" begin
        @testset "Vallado 2013 - Elliptical orbit" begin
            # Example 2-7: standard elliptical orbit
            a = 8788.1  # km
            e = 0.1712
            i = 87.87 * π/180  # rad
            Ω = 227.89 * π/180  # rad
            ω = 53.38 * π/180  # rad
            M = 92.335 * π/180  # rad (mean anomaly)
            
            # Convert to J2 equinoctial elements
            IOE = koe2IOE([a, e, i, Ω, ω, M])
            
            # Verify semi-major axis preserved
            @test IOE[1] ≈ a atol=1e-10
            
            # Round-trip: IOE -> koe should recover original
            koe_recovered = IOE2koe(IOE)
            @test koe_recovered[1] ≈ a atol=1e-10 rtol=1e-12
            @test koe_recovered[2] ≈ e atol=1e-10 rtol=1e-12
            @test koe_recovered[3] ≈ i atol=1e-10 rtol=1e-12
            @test koe_recovered[4] ≈ Ω atol=1e-10 rtol=1e-12
            @test koe_recovered[5] ≈ ω atol=1e-10 rtol=1e-12
            @test koe_recovered[6] ≈ M atol=1e-10 rtol=1e-12
        end
        
        @testset "Circular orbit (e=0)" begin
            # Edge case: circular orbit (argument of periapsis undefined)
            a = 7000.0  # km
            e = 0.0
            i = 45.0 * π/180
            Ω = 30.0 * π/180
            ω = 0.0  # Undefined for circular, but set to 0
            M = 90.0 * π/180
            
            IOE = koe2IOE([a, e, i, Ω, ω, M])
            koe_recovered = IOE2koe(IOE)
            
            @test koe_recovered[1] ≈ a atol=1e-10 rtol=1e-12
            @test koe_recovered[2] ≈ e atol=1e-10  # Should be exactly 0
            @test koe_recovered[3] ≈ i atol=1e-10 rtol=1e-12
        end
        
        @testset "Equatorial orbit (i=0)" begin
            # Edge case: equatorial orbit (RAAN undefined)
            a = 42164.0  # km (GEO altitude)
            e = 0.01
            i = 0.0
            Ω = 0.0  # Undefined for equatorial
            ω = 45.0 * π/180
            M = 180.0 * π/180
            
            IOE = koe2IOE([a, e, i, Ω, ω, M])
            koe_recovered = IOE2koe(IOE)
            
            @test koe_recovered[1] ≈ a atol=1e-10 rtol=1e-12
            @test koe_recovered[2] ≈ e atol=1e-10 rtol=1e-12
            @test koe_recovered[3] ≈ i atol=1e-10  # Should be exactly 0
        end
        
        @testset "Circular equatorial orbit (e=0, i=0)" begin
            # Double edge case: both undefined classical angles
            a = 7000.0
            e = 0.0
            i = 0.0
            Ω = 0.0
            ω = 0.0
            M = 45.0 * π/180
            
            IOE = koe2IOE([a, e, i, Ω, ω, M])
            koe_recovered = IOE2koe(IOE)
            
            @test koe_recovered[1] ≈ a atol=1e-10 rtol=1e-12
            @test koe_recovered[2] ≈ e atol=1e-10
            @test koe_recovered[3] ≈ i atol=1e-10
        end
        
        @testset "Highly eccentric orbit (e→1)" begin
            # Near-parabolic orbit tests numerical stability
            a = 20000.0
            e = 0.95
            i = 60.0 * π/180
            Ω = 120.0 * π/180
            ω = 270.0 * π/180
            M = 30.0 * π/180
            
            IOE = koe2IOE([a, e, i, Ω, ω, M])
            koe_recovered = IOE2koe(IOE)
            
            @test koe_recovered[1] ≈ a atol=1e-8 rtol=1e-10
            @test koe_recovered[2] ≈ e atol=1e-10 rtol=1e-10
            @test koe_recovered[3] ≈ i atol=1e-10 rtol=1e-12
        end
        
        @testset "Retrograde orbit (i>90°)" begin
            a = 8000.0
            e = 0.2
            i = 120.0 * π/180  # Retrograde
            Ω = 45.0 * π/180
            ω = 90.0 * π/180
            M = 180.0 * π/180
            
            IOE = koe2IOE([a, e, i, Ω, ω, M])
            koe_recovered = IOE2koe(IOE)
            
            @test koe_recovered[1] ≈ a atol=1e-10 rtol=1e-12
            @test koe_recovered[2] ≈ e atol=1e-10 rtol=1e-12
            @test koe_recovered[3] ≈ i atol=1e-10 rtol=1e-12
        end
    end
    
    @testset "State vector to J2EqOE conversions" begin
        @testset "Round-trip: rv → IOE → rv" begin
            # Test that state vector conversion preserves orbital state
            r = [7000.0, 2000.0, 1000.0]  # km
            v = [0.5, 7.0, 2.0]  # km/s
            
            # Convert to J2 equinoctial
            IOE = rv2IOE(r, v, μ)
            
            # Convert back to state vectors
            r_recovered, v_recovered = IOE2rv(IOE, μ)
            
            # Should recover original state within numerical precision
            @test r_recovered ≈ r atol=1e-8 rtol=1e-10
            @test v_recovered ≈ v atol=1e-8 rtol=1e-10
        end
        
        @testset "Circular orbit state vectors" begin
            # Circular orbit at 7000 km altitude
            r_mag = 7000.0
            v_circ = sqrt(μ / r_mag)  # Circular velocity
            r = [r_mag, 0.0, 0.0]
            v = [0.0, v_circ, 0.0]
            
            IOE = rv2IOE(r, v, μ)
            
            # Check eccentricity is near zero
            koe = IOE2koe(IOE)
            @test koe[2] ≈ 0.0 atol=1e-10
            
            # Round-trip
            r_rec, v_rec = IOE2rv(IOE, μ)
            @test r_rec ≈ r atol=1e-8 rtol=1e-10
            @test v_rec ≈ v atol=1e-8 rtol=1e-10
        end
        
        @testset "Highly elliptical orbit" begin
            # GTO-like orbit
            r_perigee = 6678.0  # km (300 km altitude)
            r_apogee = 42164.0  # km (GEO altitude)
            a = (r_perigee + r_apogee) / 2
            e = (r_apogee - r_perigee) / (r_apogee + r_perigee)
            
            # Position at perigee
            r = [r_perigee, 0.0, 0.0]
            v_perigee = sqrt(μ * (2/r_perigee - 1/a))
            v = [0.0, v_perigee, 0.0]
            
            IOE = rv2IOE(r, v, μ)
            r_rec, v_rec = IOE2rv(IOE, μ)
            
            @test r_rec ≈ r atol=1e-6 rtol=1e-10
            @test v_rec ≈ v atol=1e-6 rtol=1e-10
        end
    end
    
    @testset "Mean anomaly conversions" begin
        @testset "M ↔ ν across full range [0, 2π]" begin
            e = 0.3
            test_M = [0.0, π/6, π/4, π/3, π/2, 2π/3, π, 4π/3, 3π/2, 5π/3, 2π-0.01]
            
            for M in test_M
                # M → true anomaly → M
                # This tests KeplerSolver and the inverse transformation
                a = 10000.0
                i = 45.0 * π/180
                Ω = 30.0 * π/180
                ω = 60.0 * π/180
                
                koe1 = [a, e, i, Ω, ω, M]
                IOE = koe2IOE(koe1)
                koe2 = IOE2koe(IOE)
                
                # Mean anomaly should be preserved (modulo 2π)
                M_diff = mod(koe2[6] - M, 2π)
                @test M_diff < 1e-10 || M_diff > 2π - 1e-10
            end
        end
        
        @testset "KeplerSolver convergence" begin
            # Test Kepler equation solver for various eccentricities
            test_e = [0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99]
            M_test = [0.5, π/2, π, 3π/2]
            
            for e in test_e
                for M in M_test
                    a = 10000.0
                    i = 30.0 * π/180
                    Ω = 45.0 * π/180
                    ω = 90.0 * π/180
                    
                    koe = [a, e, i, Ω, ω, M]
                    IOE = koe2IOE(koe)
                    koe_rec = IOE2koe(IOE)
                    
                    # Should converge for all eccentricities < 1
                    @test abs(koe_rec[6] - M) < 1e-8 || abs(koe_rec[6] - M - 2π) < 1e-8
                end
            end
        end
        
        @testset "Edge cases: M=0, M=π" begin
            e = 0.5
            a = 8000.0
            i = 45.0 * π/180
            Ω = 0.0
            ω = 0.0
            
            # M = 0 (periapsis)
            koe1 = [a, e, i, Ω, ω, 0.0]
            IOE1 = koe2IOE(koe1)
            koe1_rec = IOE2koe(IOE1)
            @test koe1_rec[6] ≈ 0.0 atol=1e-10
            
            # M = π (apoapsis)
            koe2 = [a, e, i, Ω, ω, π]
            IOE2 = koe2IOE(koe2)
            koe2_rec = IOE2koe(IOE2)
            @test koe2_rec[6] ≈ π atol=1e-10
        end
    end
end
