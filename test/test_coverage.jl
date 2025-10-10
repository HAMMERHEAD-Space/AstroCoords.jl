@testset "Additional Coverage Tests" begin
    μ = 3.986004415e5

    # Test J2EqOE coordinate set
    @testset "J2EqOE Coordinate Set" begin
        state = [
            -1076.225324679696,
            -6765.896364327722,
            -332.3087833503755,
            9.356857417032581,
            -3.3123476319597557,
            -1.1880157328553503,
        ]
        cart = Cartesian(state)
        j2eqoe = J2EqOE(cart, μ)

        # Test round trip
        cart_back = Cartesian(j2eqoe, μ)
        @test params(cart) ≈ params(cart_back) rtol = 1e-12

        # Test direct construction
        j2eqoe_params = params(j2eqoe)
        j2eqoe2 = J2EqOE(j2eqoe_params)
        @test params(j2eqoe) == params(j2eqoe2)

        # Test conversion to other coordinates
        kep = Keplerian(j2eqoe, μ)
        @test kep isa Keplerian

        # Test quantities
        @test orbitalNRG(j2eqoe, μ) ≈ orbitalNRG(cart, μ) rtol = 1e-12
        @test angularMomentumQuantity(j2eqoe, μ) ≈ angularMomentumQuantity(cart, μ) rtol = 1e-12
    end

    # Test edge cases for anomaly conversions
    @testset "Anomaly Edge Cases" begin
        # Circular orbit (e=0)
        @testset "Circular Orbit" begin
            M = π / 4
            e = 0.0
            E = meanAnomaly2EccentricAnomaly(M, e)
            f = meanAnomaly2TrueAnomaly(M, e)
            @test E ≈ M atol = 1e-14
            @test f ≈ M atol = 1e-14
        end

        # Near-parabolic orbit
        @testset "Near-Parabolic Orbit" begin
            M = 0.1
            e = 0.9999
            E = meanAnomaly2EccentricAnomaly(M, e)
            f = meanAnomaly2TrueAnomaly(M, e)
            @test !isnan(E)
            @test !isnan(f)
        end

        # M = 0
        @testset "M = 0" begin
            e = 0.5
            E = meanAnomaly2EccentricAnomaly(0.0, e)
            f = meanAnomaly2TrueAnomaly(0.0, e)
            @test E ≈ 0.0 atol = 1e-14
            @test f ≈ 0.0 atol = 1e-14
        end

        # M = 2π
        @testset "M = 2π" begin
            e = 0.5
            E = meanAnomaly2EccentricAnomaly(2π, e)
            f = meanAnomaly2TrueAnomaly(2π, e)
            @test E ≈ 2π atol = 1e-12
            @test abs(f - 2π) < 1e-12 || abs(f) < 1e-12  # Could wrap to 0
        end
    end

    # Test utility functions
    @testset "Utility Function Edge Cases" begin
        @testset "angle_between_vectors" begin
            # Parallel vectors
            v1 = [1.0, 0.0, 0.0]
            v2 = [2.0, 0.0, 0.0]
            angle = angle_between_vectors(v1, v2)
            @test angle ≈ 0.0 atol = 1e-14

            # Antiparallel vectors
            v3 = [-1.0, 0.0, 0.0]
            angle2 = angle_between_vectors(v1, v3)
            @test angle2 ≈ π atol = 1e-14

            # Orthogonal vectors
            v4 = [0.0, 1.0, 0.0]
            angle3 = angle_between_vectors(v1, v4)
            @test angle3 ≈ π / 2 atol = 1e-14

            # General case
            v5 = [1.0, 1.0, 0.0]
            angle4 = angle_between_vectors(v1, v5)
            @test angle4 ≈ π / 4 atol = 1e-14
        end
    end

    # Test attitude transformation edge cases
    @testset "Attitude Transformation Edge Cases" begin
        @testset "EP2MRP Edge Cases" begin
            # Standard case
            β1 = [cos(π / 8), sin(π / 8), 0.0, 0.0]
            σ1 = EP2MRP(β1)
            β1_back = MRP2EP(σ1)
            @test β1 ≈ β1_back atol = 1e-14

            # β0 close to -1 (but not exactly -1 to avoid singularity)
            β2 = [-0.999, 0.04472135954999579, 0.0, 0.0]
            β2 = β2 / norm(β2)  # Normalize
            σ2 = EP2MRP(β2)
            @test !any(isnan.(σ2))
            @test !any(isinf.(σ2))

            # β0 = 1 (identity rotation)
            β3 = [1.0, 0.0, 0.0, 0.0]
            σ3 = EP2MRP(β3)
            @test all(abs.(σ3) .< 1e-14)
        end

        @testset "MRP2EP Edge Cases" begin
            # Small MRP
            σ1 = [0.1, 0.0, 0.0]
            β1 = MRP2EP(σ1)
            σ1_back = EP2MRP(β1)
            @test σ1 ≈ σ1_back atol = 1e-14

            # Large MRP
            σ2 = [5.0, 0.0, 0.0]
            β2 = MRP2EP(σ2)
            @test norm(β2) ≈ 1.0 atol = 1e-14

            # Zero MRP (identity)
            σ3 = [0.0, 0.0, 0.0]
            β3 = MRP2EP(σ3)
            @test β3[1] ≈ 1.0 atol = 1e-14
            @test all(abs.(β3[2:4]) .< 1e-14)
        end
    end

    # Test more coordinate pair transformations
    @testset "Additional Coordinate Pair Transformations" begin
        state = [
            -1076.225324679696,
            -6765.896364327722,
            -332.3087833503755,
            9.356857417032581,
            -3.3123476319597557,
            -1.1880157328553503,
        ]
        cart = Cartesian(state)

        pairs_to_test = [
            (Delaunay, Milankovich),
            (Keplerian, ModEq),
            (Spherical, Cylindrical),
            (USM6, USM7),
            (USMEM, Keplerian),
            (J2EqOE, Delaunay),
        ]

        for (Type1, Type2) in pairs_to_test
            state1 = Type1(cart, μ)
            state2 = Type2(state1, μ)
            state1_back = Type1(state2, μ)
            @test params(state1) ≈ params(state1_back) rtol = 1e-12
        end
    end

    # Test regularized coordinate edge cases
    @testset "Regularized Coordinate Edge Cases" begin
        state = [
            -1076.225324679696,
            -6765.896364327722,
            -332.3087833503755,
            9.356857417032581,
            -3.3123476319597557,
            -1.1880157328553503,
        ]
        cart = Cartesian(state)

        # Test with zero W
        @testset "Zero W" begin
            config = RegularizedCoordinateConfig(state, μ; W=0.0)
            @test config.W == 0.0
        end

        # Test with non-zero W
        @testset "Non-zero W" begin
            config = RegularizedCoordinateConfig(state, μ; W=1e-8)
            @test config.W == 1e-8
        end

        # Test compute_characteristic_scales
        @testset "Characteristic Scales" begin
            DU, TU = compute_characteristic_scales(state, μ)
            @test DU > 0
            @test TU > 0
            @test TU ≈ sqrt(DU^3 / μ) rtol = 1e-14
        end

        # Test compute_initial_phi
        @testset "Initial Phi" begin
            config = RegularizedCoordinateConfig(state, μ)
            ϕ = compute_initial_phi(state, μ, config)
            @test !isnan(ϕ)
            @test isfinite(ϕ)
        end
    end
end
