@testset "Test Anomaly Conversions" begin
    @testset "Basic Round-Trips" begin
        M = 2π * rand()
        e = rand()

        E = meanAnomaly2EccentricAnomaly(M, e)
        f = meanAnomaly2TrueAnomaly(M, e)

        E2 = trueAnomaly2EccentricAnomaly(f, e)
        f2 = eccentricAnomaly2TrueAnomaly(E, e)

        @test E ≈ E2 atol = 1e-14
        @test f ≈ f2 atol = 1e-14

        M_from_E = eccentricAnomaly2MeanAnomaly(E, e)
        M_from_f = trueAnomaly2MeanAnomaly(f, e)

        @test M_from_E ≈ M atol = 1e-14
        @test M_from_f ≈ M atol = 1e-14
    end

    @testset "KeplerSolver for Elliptic Orbits (e<1)" begin
        # Test various mean anomalies and eccentricities
        test_cases = [
            (M=0.0, e=0.0, desc="Circular orbit at M=0"),
            (M=π/2, e=0.0, desc="Circular orbit at M=π/2"),
            (M=0.0, e=0.1, desc="Low eccentricity at periapsis"),
            (M=1π, e=0.1, desc="Low eccentricity at apoapsis"),
            (M=π/4, e=0.5, desc="Moderate eccentricity"),
            (M=π/2, e=0.8, desc="High eccentricity"),
            (M=3π/2, e=0.5, desc="Third quadrant"),
        ]

        for tc in test_cases
            f = KeplerSolver(tc.M, tc.e)
            @test f isa Float64
            @test 0.0 ≤ f ≤ 2π  # True anomaly should be in [0, 2π]

            # Verify round-trip: f → E → M should return original M
            E = trueAnomaly2EccentricAnomaly(f, tc.e)
            M_recovered = eccentricAnomaly2MeanAnomaly(E, tc.e)
            @test M_recovered ≈ tc.M atol=1e-14
        end
    end

    @testset "Newton-Raphson iteration convergence" begin
        # Test that solver converges for various tolerances
        M = π/3
        e = 0.5

        # Default tolerance
        f1 = KeplerSolver(M, e)
        @test f1 isa Float64

        # Tighter tolerance
        f2 = KeplerSolver(M, e; tol=1e-12)
        @test f2 ≈ f1 atol=1e-11

        # Looser tolerance
        f3 = KeplerSolver(M, e; tol=1e-6)
        @test f3 ≈ f1 atol=1e-5
    end

    @testset "KeplerSolver for Hyperbolic Orbits (e>1)" begin
        # Test hyperbolic cases
        test_cases = [
            (M=0.0, e=1.1, desc="Slightly hyperbolic at M=0"),
            (M=0.5, e=1.5, desc="Moderately hyperbolic"),
            (M=1.0, e=2.0, desc="Highly hyperbolic"),
            (M=2.5, e=1.5, desc="Negative mean anomaly"),
        ]

        for tc in test_cases
            f = KeplerSolver(tc.M, tc.e)
            @test f isa Float64
            @test 0.0 ≤ f ≤ 2π

            # Verify round-trip for hyperbolic orbit
            F = trueAnomaly2EccentricAnomaly(f, tc.e)  # Hyperbolic eccentric anomaly
            M_recovered = eccentricAnomaly2MeanAnomaly(F, tc.e)
            @test M_recovered ≈ tc.M atol=1e-14
        end
    end

    @testset "Near-Parabolic Orbits (e≈0.999, e≈1.001)" begin
        # Test near-parabolic elliptic
        M = π/4
        e_near_para_ell = 0.9999
        f_ell = KeplerSolver(M, e_near_para_ell)
        @test f_ell isa Float64
        @test 0.0 ≤ f_ell ≤ 2π

        # Test near-parabolic hyperbolic
        e_near_para_hyp = 1.0001
        f_hyp = KeplerSolver(M, e_near_para_hyp)
        @test f_hyp isa Float64
        @test 0.0 ≤ f_hyp ≤ 2π

        # The results should be close but not identical
        @test abs(f_ell - f_hyp) < 0.1  # Should be close but different
    end

    @testset "Eccentric → True Anomaly Conversions" begin
        # Test eccentricAnomaly2TrueAnomaly for various cases
        test_cases = [
            (E=0.0, e=0.1), (E=π/4, e=0.3), (E=π/2, e=0.5), (E=π, e=0.7), (E=3π/2, e=0.2)
        ]

        for tc in test_cases
            f = eccentricAnomaly2TrueAnomaly(tc.E, tc.e)
            @test f isa Float64

            # Verify inverse
            E_recovered = trueAnomaly2EccentricAnomaly(f, tc.e)
            @test E_recovered ≈ tc.E atol=1e-14
        end

        # Test hyperbolic case
        F = 1.0  # Hyperbolic eccentric anomaly
        e = 1.5
        f_hyp = eccentricAnomaly2TrueAnomaly(F, e)
        @test f_hyp isa Float64
    end

    @testset "True → Eccentric → Mean Round-Trips" begin
        # Test full round-trip conversions
        test_cases = [
            (f=0.0, e=0.1),
            (f=π/6, e=0.3),
            (f=π/3, e=0.5),
            (f=π/2, e=0.6),
            (f=2π/3, e=0.4),
            (f=1π, e=0.7),
            (f=4π/3, e=0.2),
            (f=5π/3, e=0.5),
        ]

        for tc in test_cases
            # f → E → M → f
            E = trueAnomaly2EccentricAnomaly(tc.f, tc.e)
            M = eccentricAnomaly2MeanAnomaly(E, tc.e)
            f_recovered = meanAnomaly2TrueAnomaly(M, tc.e)

            @test f_recovered ≈ tc.f atol=1e-14
        end
    end

    @testset "Mean → Eccentric → True Chains" begin
        # Test conversion chains starting from mean anomaly
        test_cases = [
            (M=0.0, e=0.1),
            (M=π/4, e=0.3),
            (M=π/2, e=0.5),
            (M=3π/4, e=0.6),
            (M=1π, e=0.4),
            (M=5π/4, e=0.7),
            (M=3π/2, e=0.2),
            (M=7π/4, e=0.5),
        ]

        for tc in test_cases
            # M → E
            E = meanAnomaly2EccentricAnomaly(tc.M, tc.e)
            @test E isa Float64

            # E → f
            f = eccentricAnomaly2TrueAnomaly(E, tc.e)
            @test f isa Float64

            # Verify M → f direct
            f_direct = meanAnomaly2TrueAnomaly(tc.M, tc.e)
            @test f_direct ≈ f atol=1e-14

            # Verify full round-trip M → E → f → M
            M_recovered = trueAnomaly2MeanAnomaly(f, tc.e)
            @test M_recovered ≈ tc.M atol=1e-14
        end
    end

    @testset "Various M Values (0 to 2π)" begin
        e = 0.4
        M_values = range(0.0, 2π, length=20)

        for M in M_values
            f = KeplerSolver(M, e)
            @test 0.0 ≤ f ≤ 2π

            # Verify monotonicity: as M increases, f should generally increase
            # (though not strictly monotonic due to wrapping)
            M_recovered = trueAnomaly2MeanAnomaly(f, e)
            @test M_recovered ≈ M atol=1e-14
        end
    end

    @testset "Tolerance Parameter Effects" begin
        M = π/3
        e = 0.5

        # Test with different tolerances
        tolerances = [1e-6, 1e-8, 1e-10, 1e-12]
        results = [KeplerSolver(M, e; tol=tol) for tol in tolerances]

        # All results should be close to each other
        for i in 2:length(results)
            @test results[i] ≈ results[1] atol=1e-5
        end

        # Higher precision should give more accurate round-trip
        for tol in tolerances
            f = KeplerSolver(M, e; tol=tol)
            M_recovered = trueAnomaly2MeanAnomaly(f, e)
            @test M_recovered ≈ M atol=10*tol
        end
    end

    @testset "Edge Cases and Special Values" begin
        # Test at periapsis (f=0)
        e = 0.5
        f = 0.0
        E = trueAnomaly2EccentricAnomaly(f, e)
        @test E ≈ 0.0 atol=1e-14
        M = eccentricAnomaly2MeanAnomaly(E, e)
        @test M ≈ 0.0 atol=1e-14

        # Test at apoapsis (f=π)
        f = π
        E = trueAnomaly2EccentricAnomaly(f, e)
        @test E ≈ π atol=1e-14
        M = eccentricAnomaly2MeanAnomaly(E, e)
        @test M ≈ π atol=1e-14

        # Test with e=0 (circular orbit)
        e = 0.0
        M = π/4
        f = KeplerSolver(M, e)
        @test f ≈ M atol=1e-14  # For circular orbit, f = M = E
    end
end
