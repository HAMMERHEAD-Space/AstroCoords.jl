@testset "GEqOE Transformations" begin
    @testset "GEqOE Round Trip" begin
        μ = 3.986004415e5

        # W=0 for pure Keplerian, W=1e-3 for a small perturbation
        configs = [
            RegularizedCoordinateConfig(; W=0.0), RegularizedCoordinateConfig(; W=1e-3)
        ]

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

        # Retrograde Orbit
        r_retro = 9000.0
        v_retro = sqrt(μ / r_retro)
        i_retro = 2π / 3  # 120 degrees

        x_retro = r_retro
        vx_retro = 0.0
        vy_retro = v_retro * cos(i_retro)
        vz_retro = v_retro * sin(i_retro)
        state_retrograde = [x_retro, 0.0, 0.0, vx_retro, vy_retro, vz_retro]

        # GEqOE requires negative total energy (bound orbits only)
        states = [
            state_elliptical,
            state_circular_equatorial,
            state_ecc_equatorial,
            state_circ_inclined,
            state_retrograde,
        ]

        for config in configs
            for state in states
                base_cart_state = Cartesian(state)

                @testset "Roundtrip from Cartesian" begin
                    geqoe_state = GEqOE(base_cart_state, μ, config)
                    roundtrip_state = Cartesian(geqoe_state, μ, config)

                    @test params(roundtrip_state) ≈ params(base_cart_state) atol = 1e-8
                end
            end
        end
    end

    @testset "Element properties" begin
        μ = 3.986004415e5
        config = RegularizedCoordinateConfig(; W=0.0)

        state = [9000.0, 4000.0, 3000.0, -2.0, 4.0, 3.0]
        geqoe_vec = AstroCoords.cart2geqoe(state, μ, config)
        geqoe = GEqOE(geqoe_vec)

        # Verify all 6 elements are finite
        for i in 1:6
            @test isfinite(geqoe_vec[i])
        end

        # Computed properties
        @test geqoe.g ≈ sqrt(geqoe_vec[2]^2 + geqoe_vec[3]^2)
        @test geqoe.Ψ ≈ atan(geqoe_vec[2], geqoe_vec[3])

        # Derived quantities
        a = AstroCoords.computed_a(geqoe, μ)
        @test a > 0
        ρ = AstroCoords.computed_ρ(geqoe, μ)
        @test ρ > 0
    end

    @testset "L₀ derived quantity" begin
        μ = 3.986004415e5
        config = RegularizedCoordinateConfig(; W=0.0)
        t_test = 50.0

        state = [7000.0, 0.0, 0.0, 0.0, sqrt(μ / 7000.0), 0.0]
        geqoe_vec = AstroCoords.cart2geqoe(state, μ, config)
        geqoe = GEqOE(geqoe_vec)

        @test AstroCoords.L₀(geqoe, t_test) ≈ geqoe.L - geqoe.ν * t_test
        @test AstroCoords.L₀(geqoe, 0.0) ≈ geqoe.L
    end
end
