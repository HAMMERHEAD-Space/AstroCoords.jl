using Test
using AstroCoords

#@lit {citation="Vallado 2013", id="vallado2013", ref="Algorithm 2: Kepler's Equation"}
@testset "Anomaly Conversions Coverage" begin
    @testset "KeplerSolver - elliptic orbits (e < 1)" begin
        # Test with known values from Vallado 2013, Example 2-1
        e_ellip = 0.4
        
        # M = 0 → f = 0
        f = KeplerSolver(0.0, e_ellip)
        @test f ≈ 0.0 atol=1e-10
        
        # M = π → f = π (apoapsis)
        f_apo = KeplerSolver(π, e_ellip)
        @test f_apo ≈ π atol=1e-10
        
        # M = π/2 → intermediate value
        M_half = π/2
        f_half = KeplerSolver(M_half, e_ellip)
        @test 0.0 < f_half < π
        
        # Small eccentricity (nearly circular)
        e_small = 0.01
        M_test = π/3
        f_small = KeplerSolver(M_test, e_small)
        @test f_small ≈ M_test atol=0.02  # Should be close to M for small e
        
        # Convergence test: custom tolerance
        f_tight = KeplerSolver(M_test, e_ellip, tol=1e-14)
        @test f_tight isa Float64
        
        # Type stability
        f_float32 = KeplerSolver(Float32(M_test), e_ellip)
        @test f_float32 isa Float32
    end
    
    @testset "KeplerSolver - near-parabolic orbits (e ≈ 1)" begin
        # Just below parabolic
        e_near_para = 0.9999
        M_test = 0.1
        f = KeplerSolver(M_test, e_near_para)
        @test f > M_test  # True anomaly grows faster for high e
        @test !isnan(f)
        @test !isinf(f)
        
        # Just above parabolic
        e_near_hyp = 1.0001
        f_hyp = KeplerSolver(M_test, e_near_hyp)
        @test !isnan(f_hyp)
        @test !isinf(f_hyp)
    end
    
    @testset "KeplerSolver - hyperbolic orbits (e > 1)" begin
        e_hyp = 1.5
        
        # M = 0 → f = 0
        f_zero = KeplerSolver(0.0, e_hyp)
        @test f_zero ≈ 0.0 atol=1e-10
        
        # Positive mean anomaly
        M_pos = 1.0
        f_pos = KeplerSolver(M_pos, e_hyp)
        @test f_pos > 0.0
        @test f_pos < π  # Should be less than π for hyperbolic
        
        # Negative mean anomaly (should wrap to positive f)
        M_neg = -1.0
        f_neg = KeplerSolver(M_neg, e_hyp)
        @test f_neg > π  # Wrapped to [0, 2π)
        
        # Type stability
        f_mixed = KeplerSolver(Float64(M_pos), Float32(e_hyp))
        @test f_mixed isa Float64
    end
    
    @testset "trueAnomaly2MeanAnomaly - elliptic" begin
        e = 0.3
        
        # f = 0 → M = 0
        M = trueAnomaly2MeanAnomaly(0.0, e)
        @test M ≈ 0.0 atol=1e-12
        
        # f = π → M = π
        M_apo = trueAnomaly2MeanAnomaly(π, e)
        @test M_apo ≈ π atol=1e-12
        
        # Round-trip conversion
        f_test = π/3
        M_test = trueAnomaly2MeanAnomaly(f_test, e)
        f_back = KeplerSolver(M_test, e)
        @test f_back ≈ f_test atol=1e-10
    end
    
    @testset "trueAnomaly2MeanAnomaly - hyperbolic" begin
        e = 1.8
        
        # f = 0 → M = 0
        M = trueAnomaly2MeanAnomaly(0.0, e)
        @test M ≈ 0.0 atol=1e-12
        
        # Round-trip for hyperbolic
        f_test = π/4
        M_test = trueAnomaly2MeanAnomaly(f_test, e)
        f_back = KeplerSolver(M_test, e)
        @test f_back ≈ f_test atol=1e-10
    end
    
    @testset "trueAnomaly2EccentricAnomaly - elliptic" begin
        e = 0.4
        
        # f = 0 → E = 0
        E = trueAnomaly2EccentricAnomaly(0.0, e)
        @test E ≈ 0.0 atol=1e-12
        
        # f = π → E = π
        E_apo = trueAnomaly2EccentricAnomaly(π, e)
        @test E_apo ≈ π atol=1e-12
        
        # Intermediate value
        f_test = π/3
        E_test = trueAnomaly2EccentricAnomaly(f_test, e)
        @test 0.0 < E_test < π
        @test E_test > f_test  # E > f for 0 < f < π when e > 0
        
        # Round-trip via mean anomaly
        M = eccentricAnomaly2MeanAnomaly(E_test, e)
        f_back = eccentricAnomaly2TrueAnomaly(E_test, e)
        @test f_back ≈ f_test atol=1e-10
    end
    
    @testset "trueAnomaly2EccentricAnomaly - hyperbolic" begin
        e = 1.6
        
        # f = 0 → F = 0 (hyperbolic eccentric anomaly)
        F = trueAnomaly2EccentricAnomaly(0.0, e)
        @test F ≈ 0.0 atol=1e-12
        
        # Intermediate value
        f_test = π/6
        F_test = trueAnomaly2EccentricAnomaly(f_test, e)
        @test !isnan(F_test)
        
        # Round-trip
        f_back = eccentricAnomaly2TrueAnomaly(F_test, e)
        @test f_back ≈ f_test atol=1e-10
    end
    
    @testset "eccentricAnomaly2MeanAnomaly - elliptic" begin
        e = 0.35
        
        # E = 0 → M = 0
        M = eccentricAnomaly2MeanAnomaly(0.0, e)
        @test M ≈ 0.0 atol=1e-12
        
        # E = π → M = π
        M_apo = eccentricAnomaly2MeanAnomaly(π, e)
        @test M_apo ≈ π atol=1e-12
        
        # Kepler's equation: M = E - e*sin(E)
        E_test = π/4
        M_test = eccentricAnomaly2MeanAnomaly(E_test, e)
        expected_M = E_test - e * sin(E_test)
        @test M_test ≈ expected_M atol=1e-12
    end
    
    @testset "eccentricAnomaly2MeanAnomaly - hyperbolic" begin
        e = 1.7
        
        # F = 0 → M = 0
        M = eccentricAnomaly2MeanAnomaly(0.0, e)
        @test M ≈ 0.0 atol=1e-12
        
        # Hyperbolic Kepler equation: M = e*sinh(F) - F
        F_test = 0.5
        M_test = eccentricAnomaly2MeanAnomaly(F_test, e)
        expected_M = e * sinh(F_test) - F_test
        @test M_test ≈ expected_M atol=1e-12
    end
    
    @testset "eccentricAnomaly2TrueAnomaly - elliptic" begin
        e = 0.25
        
        # E = 0 → f = 0
        f = eccentricAnomaly2TrueAnomaly(0.0, e)
        @test f ≈ 0.0 atol=1e-12
        
        # E = π → f = π
        f_apo = eccentricAnomaly2TrueAnomaly(π, e)
        @test f_apo ≈ π atol=1e-12
        
        # Round-trip
        E_test = π/3
        f_test = eccentricAnomaly2TrueAnomaly(E_test, e)
        E_back = trueAnomaly2EccentricAnomaly(f_test, e)
        @test E_back ≈ E_test atol=1e-10
    end
    
    @testset "eccentricAnomaly2TrueAnomaly - hyperbolic" begin
        e = 1.4
        
        # F = 0 → f = 0
        f = eccentricAnomaly2TrueAnomaly(0.0, e)
        @test f ≈ 0.0 atol=1e-12
        
        # Round-trip
        F_test = 0.8
        f_test = eccentricAnomaly2TrueAnomaly(F_test, e)
        F_back = trueAnomaly2EccentricAnomaly(f_test, e)
        @test F_back ≈ F_test atol=1e-10
    end
    
    @testset "meanAnomaly2TrueAnomaly - wrapper" begin
        e = 0.3
        M_test = π/4
        
        # Should give same result as KeplerSolver
        f1 = meanAnomaly2TrueAnomaly(M_test, e)
        f2 = KeplerSolver(M_test, e)
        @test f1 ≈ f2 atol=1e-14
        
        # Custom tolerance
        f_tight = meanAnomaly2TrueAnomaly(M_test, e, tol=1e-15)
        @test f_tight ≈ f1 atol=1e-13
    end
    
    @testset "meanAnomaly2EccentricAnomaly - wrapper" begin
        e = 0.45
        M_test = π/6
        
        # Should be consistent with other conversions
        E = meanAnomaly2EccentricAnomaly(M_test, e)
        M_back = eccentricAnomaly2MeanAnomaly(E, e)
        @test M_back ≈ M_test atol=1e-10
        
        # Custom tolerance
        E_tight = meanAnomaly2EccentricAnomaly(M_test, e, tol=1e-15)
        @test E_tight ≈ E atol=1e-13
    end
    
    @testset "Full round-trip conversions" begin
        e_ellip = 0.3
        f_orig = 2.1  # rad
        
        # f → M → f
        M = trueAnomaly2MeanAnomaly(f_orig, e_ellip)
        f_back = meanAnomaly2TrueAnomaly(M, e_ellip)
        @test f_back ≈ f_orig atol=1e-10
        
        # f → E → f
        E = trueAnomaly2EccentricAnomaly(f_orig, e_ellip)
        f_back2 = eccentricAnomaly2TrueAnomaly(E, e_ellip)
        @test f_back2 ≈ f_orig atol=1e-10
        
        # E → M → E
        M2 = eccentricAnomaly2MeanAnomaly(E, e_ellip)
        E_back = meanAnomaly2EccentricAnomaly(M2, e_ellip)
        @test E_back ≈ E atol=1e-10
    end
end
