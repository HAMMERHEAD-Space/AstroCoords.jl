@testset "Coordinate Changes" begin
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

    # states = [state_elliptical,state_circular_equatorial, state_ecc_equatorial, state_ecc_equatorial, state_near_parabolic, state_hyperbolic, state_retrograde]
    states = [state_elliptical]

    @testset "Round Trip Coordinate Changes" begin
        for state in states
            cart_state = Cartesian(state)

            for coord in filter(
                T -> T ∉ (EDromo, KustaanheimoStiefel, StiefelScheifele),
                AstroCoords.COORD_TYPES,
            )
                coord_state = coord(cart_state, μ)
                for coord2 in filter(
                    T -> T ∉ (EDromo, KustaanheimoStiefel, StiefelScheifele),
                    AstroCoords.COORD_TYPES,
                )
                    coord_state2 = coord2(coord_state, μ)
                    coord_state_round_trip2 = coord(coord_state2, μ)
                    @test params(coord_state) ≈ params(coord_state_round_trip2) rtol = 1e-12
                end
            end
        end
    end
end
