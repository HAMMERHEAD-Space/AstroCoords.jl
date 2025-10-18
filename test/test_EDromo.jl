
@testset "EDromo Transformations" begin
    @testset "EDromo Round Trip" begin
        μ = 3.986004415e5

        state_elliptical = [
            -1076.225324679696
            -6765.896364327722
            -332.3087833503755
            9.356857417032581
            -3.3123476319597557
            -1.1880157328553503
        ]

        #TODO: Uncomment these as Unbounded Orbits and Singularities are resolved
        # r = 7000.0
        # v = sqrt(μ / r)

        # # Exactly equatorial (i = 0)
        # state_circular_equatorial = [r, 0.0, 0.0, 0.0, v, 0.0]

        # # Eccentric equatorial orbit
        # a_ecc = 10000.0
        # e_ecc = 0.3
        # i_ecc = 0.0
        # Ω_ecc = 1.2
        # ω_ecc = 0.5
        # f_ecc = π / 4  # true anomaly
        # koe_ecc = [a_ecc, e_ecc, i_ecc, Ω_ecc, ω_ecc, f_ecc]
        # state_ecc_equatorial = AstroCoords.koe2cart(koe_ecc, μ)

        # # Circular inclined orbit
        # a_inc = 10000.0
        # e_inc = 0.0
        # i_inc = 0.5
        # Ω_inc = 1.2
        # ω_inc = 0.5
        # f_inc = π / 4  # true anomaly
        # koe_inc = [a_inc, e_inc, i_inc, Ω_inc, ω_inc, f_inc]
        # state_ecc_equatorial = AstroCoords.koe2cart(koe_inc, μ)

        # # Near Parabolic Orbit
        # a_para = 20000.0
        # e_para = 0.99999
        # i_para = 0.4
        # Ω_para = 1.2
        # ω_para = 0.22
        # f_para = π / 6
        # koe_para = [a_para, e_para, i_para, Ω_para, ω_para, f_para]
        # state_near_parabolic = AstroCoords.koe2cart(koe_para, μ)

        # # Hyperbolic Orbit
        # a_hyp = -15000.0
        # e_hyp = 1.5
        # i_hyp = 0.5
        # Ω_hyp = 1.2
        # ω_hyp = 0.5
        # f_hyp = π / 4
        # koe_hyp = [a_hyp, e_hyp, i_hyp, Ω_hyp, ω_hyp, f_hyp]
        # state_hyperbolic = AstroCoords.koe2cart(koe_hyp, μ)

        # # Retrograde Orbit
        # r_retro = 9000.0
        # v_retro = sqrt(μ / r_retro)
        # i_retro = 2π / 3  # 120 degrees

        # x_retro = r_retro
        # vx_retro = 0.0
        # vy_retro = v_retro * cos(i_retro)
        # vz_retro = v_retro * sin(i_retro)
        # state_retrograde = [x_retro, 0.0, 0.0, vx_retro, vy_retro, vz_retro]

        # states = [state_elliptical, state_circular_equatorial, state_ecc_equatorial, state_ecc_equatorial, state_near_parabolic, state_hyperbolic, state_retrograde]
        states = [state_elliptical]

        for state in states
            base_cart_state = Cartesian(state)

            # All coordinate types except EDromo itself
            coord_types_to_test = filter(
                T -> T ∉ (EDromo, KustaanheimoStiefel, StiefelScheifele),
                AstroCoords.COORD_TYPES,
            )

            # Parameter sets to test
            time_flags = [PhysicalTime(), ConstantTime(), LinearTime()]
            W_values = [0.0, -1e-8, 1e-8]
            t0_values = [0.0, 100.0]
            ϕ_values = [0.0, 100.0]

            for FromCoord in coord_types_to_test
                @testset "Roundtrip from $FromCoord" begin
                    # Create the initial state for this coordinate type
                    from_state = FromCoord(base_cart_state, μ)

                    for flag in time_flags
                        for W in W_values
                            for t₀ in t0_values
                                for ϕ in ϕ_values
                                    @testset "Params: time_flag=$(typeof(flag)), W=$W, t₀=$t₀, ϕ=$ϕ" begin
                                        # Get the full set of parameters for this test case
                                        defaults = RegularizedCoordinateConfig(
                                            state, μ; W=W, t₀=t₀, flag_time=flag
                                        )

                                        # Perform the round trip
                                        edromo_state = EDromo(from_state, μ, ϕ, defaults)
                                        roundtrip_state = FromCoord(
                                            edromo_state, μ, ϕ, defaults
                                        )

                                        # Test for numerical equality
                                        @test params(roundtrip_state) ≈ params(from_state) atol=1e-8
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "cart2EDromo with RegularizedCoordinateConfig and ϕ" begin
        μ = 398600.4418  # Earth μ [km³/s²]
        # Test lines 20-24: cart2EDromo with ϕ parameter
        r = 8000.0  # km
        v = sqrt(μ / r)  # circular velocity
        state = [r, 0.0, 0.0, 0.0, v, 0.0]

        config = RegularizedCoordinateConfig(
            state, μ; W=0.0, t₀=0.0, flag_time=LinearTime()
        )
        ϕ = compute_initial_phi(state, μ, config)

        edromo_vec = AstroCoords.cart2EDromo(state, μ, ϕ, config)

        @test length(edromo_vec) == 8
        @test all(isfinite.(edromo_vec))
        @test edromo_vec[3] > 0  # ζ₃ = -1/(2E) > 0 for bound orbits
    end

    @testset "DU/TU Scaling" begin
        μ = 398600.4418  # Earth μ [km³/s²]
        # Test lines 26-36: non-dimensionalization
        state = [7000.0, 3000.0, 2000.0, -1.0, 5.0, 4.0]
        config = RegularizedCoordinateConfig(state, μ; W=0.0)
        ϕ = 0.5

        edromo_vec = AstroCoords.cart2EDromo(state, μ, ϕ, config)

        # Verify state is properly non-dimensionalized
        @test all(isfinite.(edromo_vec))

        # Test with different scaling
        config2 = RegularizedCoordinateConfig(
            DU=1e4, TU=100.0, W=0.0, t₀=0.0, flag_time=LinearTime()
        )
        edromo_vec2 = AstroCoords.cart2EDromo(state, μ, ϕ, config2)
        @test edromo_vec != edromo_vec2  # Different scaling
    end

    @testset "EDromo state vector construction" begin
        μ = 398600.4418  # Earth μ [km³/s²]
        # Test lines 42-62: in-plane elements computation
        r = 10000.0
        v = sqrt(μ / r)
        state = [r, 0.0, 0.0, 0.0, v, 0.0]
        config = RegularizedCoordinateConfig(state, μ)
        ϕ = compute_initial_phi(state, μ, config)

        edromo_vec = AstroCoords.cart2EDromo(state, μ, ϕ, config)

        # ζ₁, ζ₂ are in-plane elements
        @test isfinite(edromo_vec[1])  # ζ₁
        @test isfinite(edromo_vec[2])  # ζ₂
        # For circular orbit: ζ₁² + ζ₂² ≈ 0
        @test edromo_vec[1]^2 + edromo_vec[2]^2 < 1e-6

        # ζ₃ relates to energy
        @test edromo_vec[3] > 0  # ζ₃ = -1/(2E) > 0
    end

    @testset "Time transformation with flag_time" begin
        μ = 398600.4418  # Earth μ [km³/s²]
        # Test lines 67-77: different time element formulations
        state = [8000.0, 0.0, 0.0, 0.0, 7.0, 0.0]
        ϕ = 0.3

        # PhysicalTime
        config_phys = RegularizedCoordinateConfig(
            state, μ; t₀=50.0, flag_time=PhysicalTime()
        )
        edromo_phys = AstroCoords.cart2EDromo(state, μ, ϕ, config_phys)
        @test edromo_phys[8] ≈ 50.0 / config_phys.TU atol=1e-10  # ζ₈ = t₀/TU

        # ConstantTime
        config_const = RegularizedCoordinateConfig(
            state, μ; t₀=50.0, flag_time=ConstantTime()
        )
        edromo_const = AstroCoords.cart2EDromo(state, μ, ϕ, config_const)
        @test isfinite(edromo_const[8])

        # LinearTime
        config_lin = RegularizedCoordinateConfig(state, μ; t₀=50.0, flag_time=LinearTime())
        edromo_lin = AstroCoords.cart2EDromo(state, μ, ϕ, config_lin)
        @test isfinite(edromo_lin[8])
    end

    @testset "Element computation (ζ₁-ζ₈)" begin
        μ = 398600.4418  # Earth μ [km³/s²]
        # Test lines 81-117: all eight elements
        state = [9000.0, 4000.0, 3000.0, -2.0, 4.0, 3.0]
        config = RegularizedCoordinateConfig(state, μ)
        ϕ = compute_initial_phi(state, μ, config)

        edromo_vec = AstroCoords.cart2EDromo(state, μ, ϕ, config)

        # Verify all 8 elements are finite
        for i in 1:8
            @test isfinite(edromo_vec[i])
        end

        # Quaternion elements (ζ₄-ζ₇) should satisfy normalization
        q_norm = sqrt(edromo_vec[4]^2 + edromo_vec[5]^2 + edromo_vec[6]^2 + edromo_vec[7]^2)
        @test q_norm ≈ 1.0 rtol=1e-14
    end

    @testset "Inverse EDromo2cart" begin
        μ = 398600.4418  # Earth μ [km³/s²]
        # Test lines 137-189: inverse transformation
        state_orig = [7500.0, 2500.0, 1800.0, -1.5, 5.5, 3.5]
        config = RegularizedCoordinateConfig(state_orig, μ)
        ϕ = compute_initial_phi(state_orig, μ, config)

        edromo_vec = AstroCoords.cart2EDromo(state_orig, μ, ϕ, config)
        state_back = AstroCoords.EDromo2cart(edromo_vec, μ, ϕ, config)

        @test state_back[1] ≈ state_orig[1] rtol=1e-14  # x
        @test state_back[2] ≈ state_orig[2] rtol=1e-14  # y
        @test state_back[3] ≈ state_orig[3] rtol=1e-14  # z
        @test state_back[4] ≈ state_orig[4] rtol=1e-14  # vx
        @test state_back[5] ≈ state_orig[5] rtol=1e-14  # vy
        @test state_back[6] ≈ state_orig[6] rtol=1e-14  # vz
    end

    @testset "Different ϕ (Anomaly-Like) Values" begin
        μ = 398600.4418  # Earth μ [km³/s²]
        # Test transformation at different ϕ values
        state = [8000.0, 0.0, 0.0, 0.0, 7.0, 0.0]
        config = RegularizedCoordinateConfig(state, μ)

        ϕ_values = [0.0, π/4, π/2, π, 3π/2, 2π]

        for ϕ in ϕ_values
            edromo_vec = AstroCoords.cart2EDromo(state, μ, ϕ, config)
            @test all(isfinite.(edromo_vec))

            # Verify round-trip works at each ϕ
            state_back = AstroCoords.EDromo2cart(edromo_vec, μ, ϕ, config)
            @test state ≈ state_back rtol=1e-10
        end
    end
end
