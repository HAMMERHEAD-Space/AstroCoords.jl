@testset "Test Quantities" begin
    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ]

    μ = 3.986004415e5

    cart_state = Cartesian(state)

    NRG = orbitalNRG(cart_state, μ)
    h_vec = angularMomentumVector(cart_state, μ)
    h = angularMomentumQuantity(cart_state, μ)

    # Regression Tests
    @test NRG ≈ -8.146487030627135 rtol = 1e-14
    @test h_vec ≈ [6937.269116080104; -4387.938522053871; 66872.35998509153] rtol = 1e-14
    @test h ≈ 67374.26984567562 rtol = 1e-14

    coordinate_sets = [
        Cartesian,
        Delaunay,
        Keplerian,
        Milankovich,
        ModEq,
        Cylindrical,
        Spherical,
        USM7,
        USM6,
        USMEM,
    ]

    for coord in coordinate_sets
        coord_state = coord(cart_state, μ)
        NRG2 = orbitalNRG(coord_state, μ)
        h_vec2 = angularMomentumVector(coord_state, μ)
        h2 = angularMomentumQuantity(coord_state, μ)

        @test NRG ≈ NRG2 rtol = 1e-14
        @test h_vec ≈ h_vec2 rtol = 1e-14
        @test h ≈ h2 rtol = 1e-14
    end

    @testset "Test EDromo Quantities" begin
        for W in (0.0, 1e-8)
            for t₀ in (0.0, 100.0)
                for flag_time in (PhysicalTime(), ConstantTime(), LinearTime())
                    edromo_params = set_edromo_configurations(
                        state, μ; W₀=W, t₀=t₀, flag_time=flag_time
                    )

                    edromo_state = EDromo(cart_state, μ; edromo_params...)
                    NRG2 = orbitalNRG(edromo_state, μ; edromo_params...)
                    h_vec2 = angularMomentumVector(edromo_state, μ; edromo_params...)
                    h2 = angularMomentumQuantity(edromo_state, μ; edromo_params...)

                    @test NRG ≈ NRG2 rtol = 1e-12
                    @test h_vec ≈ h_vec2 rtol = 1e-12
                    @test h ≈ h2 rtol = 1e-12
                end
            end
        end
    end

    @testset "Test Kustaanheimo-Stiefel Quantities" begin
        for Vpot in (0.0, 1e-8)
            for t₀ in (0.0, 100.0)
                for flag_time in (PhysicalTime(), LinearTime())
                    ks_params = set_ks_configurations(
                        state, μ; Vpot=Vpot, t₀=t₀, flag_time=flag_time
                    )

                    ks_state = KustaanheimoStiefel(cart_state, μ; ks_params...)
                    NRG2 = orbitalNRG(ks_state, μ; ks_params...)
                    h_vec2 = angularMomentumVector(ks_state, μ; ks_params...)
                    h2 = angularMomentumQuantity(ks_state, μ; ks_params...)

                    @test NRG ≈ NRG2 rtol = 1e-12
                    @test h_vec ≈ h_vec2 rtol = 1e-12
                    @test h ≈ h2 rtol = 1e-12
                end
            end
        end
    end

    @testset "Test Stiefel-Scheifele Quantities" begin
        for W in (0.0, 1e-8)
            for t₀ in (0.0, 100.0)
                for flag_time in (PhysicalTime(), LinearTime())
                    ss_params = set_stiefelscheifele_configurations(
                        state, μ; W=W, t₀=t₀, flag_time=flag_time
                    )

                    ss_state = StiefelScheifele(cart_state, μ; ss_params...)
                    NRG2 = orbitalNRG(ss_state, μ; ss_params...)
                    h_vec2 = angularMomentumVector(ss_state, μ; ss_params...)
                    h2 = angularMomentumQuantity(ss_state, μ; ss_params...)

                    @test NRG ≈ NRG2 rtol = 1e-12
                    @test h_vec ≈ h_vec2 rtol = 1e-12
                    @test h ≈ h2 rtol = 1e-12
                end
            end
        end
    end
end
