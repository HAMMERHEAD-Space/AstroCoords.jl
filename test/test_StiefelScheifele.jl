@testset "Stiefel-Scheifele Transformations" begin
    @testset "Stiefel-Scheifele Round Trip" begin
        μ = 3.986004415e5

        state_elliptical = [
            -1076.225324679696
            -6765.896364327722
            -332.3087833503755
            9.356857417032581
            -3.3123476319597557
            -1.1880157328553503
        ]

        r = 7000.0
        v = sqrt(μ / r)

        # Exactly equatorial (i = 0)
        state_circular_equatorial = [r, 0.0, 0.0, 0.0, v, 0.0]

        # Eccentric equatorial orbit
        a_ecc = 10000.0
        e_ecc = 0.3
        i_ecc = 0.0
        Ω_ecc = 1.2
        ω_ecc = 0.5
        f_ecc = π / 4  # true anomaly
        koe_ecc = [a_ecc, e_ecc, i_ecc, Ω_ecc, ω_ecc, f_ecc]
        state_ecc_equatorial = AstroCoords.koe2cart(koe_ecc, μ)

        # Circular inclined orbit
        a_inc = 10000.0
        e_inc = 0.0
        i_inc = 0.5
        Ω_inc = 1.2
        ω_inc = 0.5
        f_inc = π / 4  # true anomaly
        koe_inc = [a_inc, e_inc, i_inc, Ω_inc, ω_inc, f_inc]
        state_circ_inclined = AstroCoords.koe2cart(koe_inc, μ)

        # Near Parabolic Orbit
        a_para = 20000.0
        e_para = 0.99999
        i_para = 0.4
        Ω_para = 1.2
        ω_para = 0.22
        f_para = π / 6
        koe_para = [a_para, e_para, i_para, Ω_para, ω_para, f_para]
        state_near_parabolic = AstroCoords.koe2cart(koe_para, μ)

        # Hyperbolic Orbit
        a_hyp = -15000.0
        e_hyp = 1.5
        i_hyp = 0.5
        Ω_hyp = 1.2
        ω_hyp = 0.5
        f_hyp = π / 4
        koe_hyp = [a_hyp, e_hyp, i_hyp, Ω_hyp, ω_hyp, f_hyp]
        state_hyperbolic = AstroCoords.koe2cart(koe_hyp, μ)

        # Retrograde Orbit
        r_retro = 9000.0
        v_retro = sqrt(μ / r_retro)
        i_retro = 2π / 3  # 120 degrees

        x_retro = r_retro
        vx_retro = 0.0
        vy_retro = v_retro * cos(i_retro)
        vz_retro = v_retro * sin(i_retro)
        state_retrograde = [x_retro, 0.0, 0.0, vx_retro, vy_retro, vz_retro]

        states = [
            state_elliptical,
            state_circular_equatorial,
            state_ecc_equatorial,
            state_circ_inclined,
            state_near_parabolic,
            state_hyperbolic,
            state_retrograde,
        ]

        for state in states
            base_cart_state = Cartesian(state)

            # Parameter sets to test
            time_flags = [PhysicalTime(), LinearTime()]
            W_values = [0.0, 1e-8]
            t0_values = [0.0, 100.0]

            @testset "Roundtrip from Cartesian" begin
                # Create the initial state for this coordinate type
                from_state = Cartesian(base_cart_state, μ)

                for flag in time_flags
                    for W in W_values
                        for t₀ in t0_values
                            @testset "Params: time_flag=$(typeof(flag)), W=$W, t₀=$t₀" begin
                                # Get the full set of parameters for this test case
                                defaults = RegularizedCoordinateConfig(
                                    state, μ; W=W, t₀=t₀, flag_time=flag
                                )

                                # Compute phi for this configuration
                                ϕ = compute_initial_phi(state, μ, defaults)

                                # Perform the round trip
                                ss_state = StiefelScheifele(base_cart_state, μ, ϕ, defaults)
                                roundtrip_state = Cartesian(ss_state, μ, ϕ, defaults)

                                # Test for numerical equality
                                @test params(roundtrip_state) ≈ params(base_cart_state) atol=1e-8
                            end
                        end
                    end
                end
            end
        end
    end

    @testset "RegularizedCoordinateConfig Usage" begin
        # Test parameters
        μ = 398600.4418  # Earth gravitational parameter (km³/s²)
        DU = 6378.137    # km
        TU = 806.8       # seconds
        W = 0.0          # no perturbing potential
        t₀ = 0.0
        # Test config creation with keyword constructor
        config1 = RegularizedCoordinateConfig(
            DU=DU, TU=TU, W=W, t₀=t₀, flag_time=PhysicalTime()
        )
        @test config1.DU ≈ DU
        @test config1.TU ≈ TU
        @test config1.W ≈ W
        @test config1.t₀ ≈ t₀
        @test config1.flag_time isa PhysicalTime

        # Test config creation from state
        cart_state = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        config2 = RegularizedCoordinateConfig(cart_state, μ, W=W, t₀=t₀)
        @test config2.DU > 0
        @test config2.TU > 0

        # Test different flag_time types
        config_const = RegularizedCoordinateConfig(
            DU=DU, TU=TU, W=W, t₀=t₀, flag_time=ConstantTime()
        )
        @test config_const.flag_time isa ConstantTime

        config_linear = RegularizedCoordinateConfig(
            DU=DU, TU=TU, W=W, t₀=t₀, flag_time=LinearTime()
        )
        @test config_linear.flag_time isa LinearTime
    end
end
