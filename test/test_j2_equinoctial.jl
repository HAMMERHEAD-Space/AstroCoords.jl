using Test
using AstroCoords
using LinearAlgebra

#@lit {citation="Vallado 2013", id="vallado2013", ref="Example 2-5, p.103"}
@testset "J2 Equinoctial Elements" begin
    # Earth gravitational parameter
    μ = 398600.4418  # km³/s²
    
    @testset "koe2IOE and IOE2koe round-trip" begin
        #@lit {citation="Vallado 2013", id="vallado2013", ref="Example 2-5, p.103"}
        @testset "Vallado 2013 — Example 2-5: Standard elliptical orbit" begin
            # Classical orbital elements from Vallado Example 2-5
            a = 8788.1      # km (semi-major axis)
            e = 0.1712      # eccentricity
            i = 87.87 * π/180  # inclination (rad)
            Ω = 227.89 * π/180 # RAAN (rad)
            ω = 53.38 * π/180  # argument of periapsis (rad)
            M = 92.335 * π/180 # mean anomaly (rad)
            
            koe_input = [a, e, i, Ω, ω, M]
            
            # Convert to IOE and back
            ioe = koe2IOE(koe_input)
            koe_output = IOE2koe(ioe)
            
            # Round-trip should preserve elements within numerical precision
            # Tolerances account for trigonometric conversions
            @test koe_output[1] ≈ koe_input[1] atol=1e-6 rtol=1e-10  # a
            @test koe_output[2] ≈ koe_input[2] atol=1e-8 rtol=1e-10  # e
            @test koe_output[3] ≈ koe_input[3] atol=1e-8 rtol=1e-10  # i
            @test koe_output[4] ≈ koe_input[4] atol=1e-8 rtol=1e-10  # Ω
            @test koe_output[5] ≈ koe_input[5] atol=1e-8 rtol=1e-10  # ω
            @test koe_output[6] ≈ koe_input[6] atol=1e-8 rtol=1e-10  # M
        end
        
        @testset "Singularity: Circular orbit (e ≈ 0)" begin
            # Near-circular orbit tests ω definition
            a = 7000.0
            e = 1e-10  # Nearly circular
            i = 45.0 * π/180
            Ω = 30.0 * π/180
            ω = 60.0 * π/180  # Should be handled despite e≈0
            M = 120.0 * π/180
            
            koe_input = [a, e, i, Ω, ω, M]
            ioe = koe2IOE(koe_input)
            koe_output = IOE2koe(ioe)
            
            # For circular orbits, ω is undefined but algorithm should handle gracefully
            @test koe_output[1] ≈ koe_input[1] atol=1e-6 rtol=1e-10
            @test koe_output[2] ≈ koe_input[2] atol=1e-12 rtol=1e-8
            @test koe_output[3] ≈ koe_input[3] atol=1e-10 rtol=1e-10
        end
        
        @testset "Singularity: Equatorial orbit (i ≈ 0)" begin
            # Near-equatorial orbit tests Ω definition
            a = 7000.0
            e = 0.1
            i = 1e-10  # Nearly equatorial
            Ω = 45.0 * π/180  # Should be handled despite i≈0
            ω = 60.0 * π/180
            M = 120.0 * π/180
            
            koe_input = [a, e, i, Ω, ω, M]
            ioe = koe2IOE(koe_input)
            koe_output = IOE2koe(ioe)
            
            # For equatorial orbits, Ω is undefined but should be handled
            @test koe_output[1] ≈ koe_input[1] atol=1e-6 rtol=1e-10
            @test koe_output[2] ≈ koe_input[2] atol=1e-10 rtol=1e-10
            @test koe_output[3] ≈ koe_input[3] atol=1e-12 rtol=1e-8
        end
        
        @testset "Singularity: Polar orbit (i ≈ π/2)" begin
            # Polar orbit
            a = 7000.0
            e = 0.05
            i = π/2  # Polar
            Ω = 45.0 * π/180
            ω = 60.0 * π/180
            M = 120.0 * π/180
            
            koe_input = [a, e, i, Ω, ω, M]
            ioe = koe2IOE(koe_input)
            koe_output = IOE2koe(ioe)
            
            @test koe_output[1] ≈ koe_input[1] atol=1e-6 rtol=1e-10
            @test koe_output[2] ≈ koe_input[2] atol=1e-10 rtol=1e-10
            @test koe_output[3] ≈ koe_input[3] atol=1e-10 rtol=1e-10
        end
    end
    
    @testset "J2EqOE propagation and conversions" begin
        @testset "IOE2J2EqOE state propagation" begin
            # Test J2 equinoctial element propagation
            a = 7000.0  # km
            e = 0.01
            i = 45.0 * π/180
            Ω = 0.0
            ω = 0.0
            M = 0.0
            
            koe = [a, e, i, Ω, ω, M]
            ioe = koe2IOE(koe)
            
            # Propagate with J2
            Re = 6378.137  # Earth radius (km)
            J2 = 1.08263e-3
            dt = 600.0  # 10 minutes
            
            j2eq = IOE2J2EqOE(ioe, μ, Re, J2, dt)
            
            # J2EqOE should be 6-element vector
            @test length(j2eq) == 6
            
            # Semi-major axis should be preserved (J2 is conservative)
            @test j2eq[1] ≈ a atol=1.0 rtol=1e-6
            
            # Eccentricity magnitude should be preserved
            h = j2eq[2]
            k = j2eq[3]
            e_prop = sqrt(h^2 + k^2)
            @test e_prop ≈ e atol=1e-6 rtol=1e-4
        end
        
        @testset "J2EqOE2IOE conversion" begin
            # Create J2 equinoctial elements
            a = 7000.0
            h = 0.01  # e*sin(ω+Ω)
            k = 0.005  # e*cos(ω+Ω)
            p = 0.3   # tan(i/2)*sin(Ω)
            q = 0.2   # tan(i/2)*cos(Ω)
            λ = 1.5   # Mean longitude
            
            j2eq = [a, h, k, p, q, λ]
            
            # Convert to IOE
            ioe = J2EqOE2IOE(j2eq)
            
            # IOE should be 6-element vector
            @test length(ioe) == 6
            
            # Semi-major axis preserved
            @test ioe[1] ≈ a atol=1e-10
            
            # Eccentricity
            e_recovered = sqrt(h^2 + k^2)
            @test abs(ioe[2]) ≈ e_recovered atol=1e-10 rtol=1e-8
        end
        
        @testset "J2EqOE2koe direct conversion" begin
            # Test direct conversion from J2EqOE to Keplerian
            a = 7000.0
            h = 0.01
            k = 0.005
            p = 0.3
            q = 0.2
            λ = 1.5
            
            j2eq = [a, h, k, p, q, λ]
            
            # Convert directly to Keplerian
            koe = J2EqOE2koe(j2eq)
            
            # Should be 6-element vector
            @test length(koe) == 6
            
            # Semi-major axis preserved
            @test koe[1] ≈ a atol=1e-10
            
            # Eccentricity from h,k
            e_expected = sqrt(h^2 + k^2)
            @test koe[2] ≈ e_expected atol=1e-10 rtol=1e-8
            
            # Inclination from p,q (should be positive)
            @test koe[3] >= 0.0
            @test koe[3] <= π
        end
        
        @testset "Critical inclination (i = 63.4°)" begin
            # At critical inclination, J2 perturbations have special behavior
            i_crit = 63.4349488 * π/180  # Critical inclination
            a = 7000.0
            e = 0.01
            Ω = 45.0 * π/180
            ω = 30.0 * π/180
            M = 0.0
            
            koe = [a, e, i_crit, Ω, ω, M]
            ioe = koe2IOE(koe)
            
            Re = 6378.137
            J2 = 1.08263e-3
            dt = 1000.0
            
            # Should handle critical inclination without issues
            j2eq = IOE2J2EqOE(ioe, μ, Re, J2, dt)
            
            @test length(j2eq) == 6
            @test isfinite(j2eq[1])  # No NaN or Inf
            @test all(isfinite.(j2eq))
        end
        
        @testset "Round-trip: koe → IOE → J2EqOE → IOE → koe" begin
            a = 7500.0
            e = 0.05
            i = 50.0 * π/180
            Ω = 20.0 * π/180
            ω = 40.0 * π/180
            M = 100.0 * π/180
            
            koe1 = [a, e, i, Ω, ω, M]
            ioe1 = koe2IOE(koe1)
            
            Re = 6378.137
            J2 = 1.08263e-3
            dt = 0.0  # No propagation, just conversion
            
            j2eq = IOE2J2EqOE(ioe1, μ, Re, J2, dt)
            ioe2 = J2EqOE2IOE(j2eq)
            koe2 = IOE2koe(ioe2)
            
            # Should preserve elements through conversions
            @test koe2[1] ≈ koe1[1] atol=1e-4 rtol=1e-8
            @test koe2[2] ≈ koe1[2] atol=1e-6 rtol=1e-6
            @test koe2[3] ≈ koe1[3] atol=1e-6 rtol=1e-6
        end
    end
end
