@testset "Kustaanheimo-Stiefel Transformations" begin
    @testset "Kustaanheimo-Stiefel Round Trip" begin
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

            coord_types_to_test = filter(
                T -> T ∉ (EDromo, KustaanheimoStiefel, StiefelScheifele),
                AstroCoords.COORD_TYPES,
            )

            # Parameter sets to test
            time_flags = [PhysicalTime(), LinearTime()]
            Vpot_values = [0.0, 1e-8]
            t0_values = [0.0, 100.0]

            for FromCoord in coord_types_to_test
                @testset "Roundtrip from $FromCoord" begin
                    # Create the initial state for this coordinate type
                    from_state = FromCoord(base_cart_state, μ)

                    for flag in time_flags
                        for Vpot in Vpot_values
                            for t₀ in t0_values
                                @testset "Params: time_flag=$(typeof(flag)), Vpot=$Vpot, t₀=$t₀" begin

                                    # Get the full set of parameters for this test case
                                    configs = RegularizedCoordinateConfig(
                                        state, μ; W=Vpot, t₀=t₀, flag_time=flag
                                    )

                                    # Perform the round trip
                                    ks_state = KustaanheimoStiefel(from_state, μ, configs)
                                    roundtrip_state = FromCoord(ks_state, μ, configs)

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

    @testset "cart2KS with RegularizedCoordinateConfig" begin
        μ = 398600.4418  # Earth μ [km³/s²]
        # Test lines 16-20: cart2KS with config
        r = 7000.0  # km
        v = sqrt(μ / r)  # circular velocity
        state = [r, 0.0, 0.0, 0.0, v, 0.0]

        config = RegularizedCoordinateConfig(
            state, μ; W=0.0, t₀=0.0, flag_time=LinearTime()
        )
        ks_vec = AstroCoords.cart2KS(state, μ, config)

        @test length(ks_vec) == 10
        @test all(isfinite.(ks_vec))
        # Verify KS position magnitude relates to physical r
        u_mag_sq = ks_vec[1]^2 + ks_vec[2]^2 + ks_vec[3]^2 + ks_vec[4]^2
        @test u_mag_sq * config.DU ≈ r rtol=1e-10
    end

    @testset "W Parameter Computation and Scaling" begin
        μ = 398600.4418  # Earth μ [km³/s²]
        # Test lines 25-33: W parameter in transformation
        r = 8000.0
        v = sqrt(μ / r)
        state = [r, 0.0, 0.0, 0.0, v, 0.0]

        # Test with zero W
        config_zero = RegularizedCoordinateConfig(state, μ; W=0.0)
        ks_zero = AstroCoords.cart2KS(state, μ, config_zero)
        @test isfinite(ks_zero[9])  # h (energy)

        # Test with non-zero W
        W_test = 1e-3
        config_W = RegularizedCoordinateConfig(state, μ; W=W_test)
        ks_W = AstroCoords.cart2KS(state, μ, config_W)
        @test isfinite(ks_W[9])
        # Energy should differ when W is non-zero
        @test abs(ks_zero[9] - ks_W[9]) > 1e-10
    end

    @testset "KS matrix transformations" begin
        μ = 398600.4418  # Earth μ [km³/s²]
        # Test lines 39-48: KS position computation
        state = [6000.0, 3000.0, 2000.0, -1.0, 5.0, 3.0]
        config = RegularizedCoordinateConfig(state, μ)
        ks_vec = AstroCoords.cart2KS(state, μ, config)

        # Verify KS coordinates follow bilinear relations
        u1, u2, u3, u4 = ks_vec[1:4]
        r_x = (u1^2 - u2^2 - u3^2 + u4^2) * config.DU
        r_y = 2 * (u1*u2 - u3*u4) * config.DU
        r_z = 2 * (u1*u3 + u2*u4) * config.DU

        @test r_x ≈ state[1] rtol=1e-10
        @test r_y ≈ state[2] rtol=1e-10
        @test r_z ≈ state[3] rtol=1e-10
    end

    @testset "Velocity transformations" begin
        μ = 398600.4418  # Earth μ [km³/s²]
        # Test lines 52-64: KS velocity computation
        state = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        config = RegularizedCoordinateConfig(state, μ)
        ks_vec = AstroCoords.cart2KS(state, μ, config)

        # KS velocity components should be finite
        @test all(isfinite.(ks_vec[5:8]))

        # Test velocity magnitude relationship
        u = SVector{4}(ks_vec[1:4])
        u_dot = SVector{4}(ks_vec[5:8])
        @test isfinite(dot(u, u_dot))
    end

    @testset "Time-element computation" begin
        μ = 398600.4418  # Earth μ [km³/s²]
        # Test lines 67-87: time element with different flag_time
        state = [8000.0, 0.0, 0.0, 0.0, 7.0, 0.0]

        # PhysicalTime
        config_phys = RegularizedCoordinateConfig(state, μ; flag_time=PhysicalTime())
        ks_phys = AstroCoords.cart2KS(state, μ, config_phys)
        @test ks_phys[10] ≈ 0.0 atol=1e-15  # τ = t₀/TU

        # LinearTime
        config_lin = RegularizedCoordinateConfig(state, μ; flag_time=LinearTime())
        ks_lin = AstroCoords.cart2KS(state, μ, config_lin)
        @test isfinite(ks_lin[10])
    end

    @testset "Inverse transformation KS2cart" begin
        μ = 398600.4418  # Earth μ [km³/s²]
        # Test lines 83-99: KS2cart function
        state_orig = [7500.0, 2000.0, 1500.0, -0.5, 6.0, 4.0]
        config = RegularizedCoordinateConfig(state_orig, μ)

        ks_vec = AstroCoords.cart2KS(state_orig, μ, config)
        state_back = AstroCoords.KS2cart(ks_vec, μ, config)

        @test state_back[1] ≈ state_orig[1] rtol=1e-10  # x
        @test state_back[2] ≈ state_orig[2] rtol=1e-10  # y
        @test state_back[3] ≈ state_orig[3] rtol=1e-10  # z
        @test state_back[4] ≈ state_orig[4] rtol=1e-10  # vx
        @test state_back[5] ≈ state_orig[5] rtol=1e-10  # vy
        @test state_back[6] ≈ state_orig[6] rtol=1e-10  # vz
    end
end
