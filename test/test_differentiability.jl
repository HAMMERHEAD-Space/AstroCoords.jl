########################################################################################
#
# Tests for Differentiation Across the Coordinate Sets
#
########################################################################################
# Currently Supported & Tested
# Diffractor, Enzyme, ForwardDiff, FiniteDiff, Mooncake, PolyesterForwardDiff, Zygote
########################################################################################
const _BACKENDS = (
    ("ForwardDiff", AutoForwardDiff()),
    ("Diffractor", AutoDiffractor()),
    ("Enzyme", AutoEnzyme()),
    ("Mooncake", AutoMooncake(; config=nothing)),
    ("PolyesterForwardDiff", AutoPolyesterForwardDiff()),
)

@testset "Coordinate Transformation Differentiation" begin
    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ]

    μ = 3.986004415e5

    for backend in _BACKENDS
        testset_name = "Coordinate Set Transformation " * backend[1]
        @testset "$testset_name" begin
            for set in filter(T -> T != EDromo, AstroCoords.COORD_TYPES)
                f_fd, df_fd = value_and_jacobian(
                    (x) -> set(Cartesian(x), μ), AutoFiniteDiff(), state
                )

                f_fd2, df_fd2 = value_and_derivative(
                    (x) -> set(Cartesian(state), x), AutoFiniteDiff(), μ
                )

                f_ad, df_ad = value_and_jacobian(
                    (x) -> Array(params(set(Cartesian(x), μ))), backend[2], state
                )

                @test f_fd == f_ad

                #TODO: Diffractor has some issue with the USM sets
                if backend[1] == "Diffractor"
                    @test df_fd ≈ df_ad rtol = 2e0
                else
                    @test df_fd ≈ df_ad atol = 1e-2
                end

                f_ad2, df_ad2 = value_and_derivative(
                    (x) -> Array(params(set(Cartesian(state), x))), backend[2], μ
                )

                @test f_fd2 == f_ad2
                @test df_fd2 ≈ df_ad2 atol = 1e-4
            end
        end
    end

    @testset "Coordinate Set Transformation Zygote" begin
        for set in filter(T -> T != EDromo, AstroCoords.COORD_TYPES)
            f_fd, df_fd = value_and_jacobian(
                (x) -> set(Cartesian(x), μ), AutoFiniteDiff(), state
            )

            f_fd2, df_fd2 = value_and_derivative(
                (x) -> set(Cartesian(state), x), AutoFiniteDiff(), μ
            )

            try
                f_ad, df_ad = value_and_jacobian(
                    (x) -> set(Cartesian(x), μ), AutoZygote(), state
                )
                @test f_fd == f_ad
                @test df_fd ≈ df_ad atol = 1e-2
            catch err
                @test err isa MethodError
                @test startswith(
                    sprint(showerror, err),
                    "MethodError: no method matching iterate(::Nothing)",
                )
            end

            try
                f_ad2, df_ad2 = value_and_derivative(
                    (x) -> set(Cartesian(state), x), AutoZygote(), μ
                )
                @test f_fd2 == f_ad2
                @test df_fd2 ≈ something.(df_ad2, 0.0) atol = 1e-4
            catch err
                @test err isa MethodError
                @test startswith(
                    sprint(showerror, err),
                    "MethodError: no method matching iterate(::Nothing)",
                )
            end
        end
    end
end

@testset "J2EqOE Step-by-Step Differentiation" begin
    μ = 3.986004415e5
    J2 = 1.0826261738522e-03
    Req = 6.378137e+03

    cart_state = Cartesian(
        -1076.225324679696,
        -6765.896364327722,
        -332.3087833503755,
        9.356857417032581,
        -3.3123476319597557,
        -1.1880157328553503,
    )
    kep_state = Array(params(Keplerian(cart_state, μ)))
    state = Array(AstroCoords.koe2IOE(kep_state, μ))

    # Test each step with Mooncake backend
    @testset "Step 1: Keplerian to J2 Elements" begin
        f_fd, df_fd = value_and_jacobian(
            (x) -> collect(AstroCoords._step1(x, μ)), AutoFiniteDiff(), state
        )

        f_ad, df_ad = value_and_jacobian(
            (x) -> collect(AstroCoords._step1(x, μ)), AutoMooncake(; config=nothing), state
        )

        @test f_fd == f_ad
        @test df_fd ≈ df_ad atol = 1e-2
    end

    @testset "Step 2: J2 Parameters" begin
        aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, fⱼ₂, Lⱼ₂ = AstroCoords._step1(state, μ)

        f_fd, df_fd = value_and_jacobian(
            (x) -> collect(AstroCoords._step2(J2, Req, x[1], x[2], x[3])),
            AutoFiniteDiff(),
            [eⱼ₂, aⱼ₂, Iⱼ₂],
        )

        f_ad, df_ad = value_and_jacobian(
            (x) -> collect(AstroCoords._step2(J2, Req, x[1], x[2], x[3])),
            AutoMooncake(; config=nothing),
            [eⱼ₂, aⱼ₂, Iⱼ₂],
        )

        @test f_fd == f_ad
        @test df_fd ≈ df_ad atol = 1e-2
    end

    @testset "Step 3: Eccentric Anomaly" begin
        aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, fⱼ₂, Lⱼ₂ = AstroCoords._step1(state, μ)

        f_fd, df_fd = value_and_jacobian(
            (x) -> collect(AstroCoords._step3(x[1], x[2])), AutoFiniteDiff(), [fⱼ₂, eⱼ₂]
        )

        f_ad, df_ad = value_and_jacobian(
            (x) -> collect(AstroCoords._step3(x[1], x[2])),
            AutoMooncake(; config=nothing),
            [fⱼ₂, eⱼ₂],
        )

        @test f_fd == f_ad
        @test df_fd ≈ df_ad atol = 1e-2
    end

    @testset "Step 4: Radius and True Anomaly" begin
        aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, fⱼ₂, Lⱼ₂ = AstroCoords._step1(state, μ)
        _, η, _, _, _ = AstroCoords._step2(J2, Req, eⱼ₂, aⱼ₂, Iⱼ₂)
        Eⱼ₂ = AstroCoords._step3(fⱼ₂, eⱼ₂)

        f_fd, df_fd = value_and_jacobian(
            (x) -> collect(AstroCoords._step4(x[1], x[2], x[3], x[4], x[5])),
            AutoFiniteDiff(),
            [aⱼ₂, eⱼ₂, η, Eⱼ₂, Lⱼ₂],
        )

        f_ad, df_ad = value_and_jacobian(
            (x) -> collect(AstroCoords._step4(x[1], x[2], x[3], x[4], x[5])),
            AutoMooncake(; config=nothing),
            [aⱼ₂, eⱼ₂, η, Eⱼ₂, Lⱼ₂],
        )

        @test f_fd == f_ad
        @test df_fd ≈ df_ad atol = 1e-2
    end

    @testset "Step 5: Semi-major Axis and RAAN Correction" begin
        aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, fⱼ₂, Lⱼ₂ = AstroCoords._step1(state, μ)
        _, η, γ, γ′, θ = AstroCoords._step2(J2, Req, eⱼ₂, aⱼ₂, Iⱼ₂)
        Eⱼ₂ = AstroCoords._step3(fⱼ₂, eⱼ₂)
        rⱼ₂, νⱼ₂, Lⱼ₂ = AstroCoords._step4(aⱼ₂, eⱼ₂, η, Eⱼ₂, Lⱼ₂)

        f_fd, df_fd = value_and_jacobian(
            (x) -> collect(
                AstroCoords._step5(
                    x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10]
                ),
            ),
            AutoFiniteDiff(),
            [aⱼ₂, γ, γ′, θ, rⱼ₂, η, eⱼ₂, gⱼ₂, νⱼ₂, Lⱼ₂],
        )

        f_ad, df_ad = value_and_jacobian(
            (x) -> collect(
                AstroCoords._step5(
                    x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10]
                ),
            ),
            AutoMooncake(; config=nothing),
            [aⱼ₂, γ, γ′, θ, rⱼ₂, η, eⱼ₂, gⱼ₂, νⱼ₂, Lⱼ₂],
        )

        @test f_fd == f_ad
        @test df_fd ≈ df_ad atol = 1e-2
    end

    @testset "Step 6: Longitude Sum" begin
        aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, fⱼ₂, Lⱼ₂ = AstroCoords._step1(state, μ)
        _, η, γ, γ′, θ = AstroCoords._step2(J2, Req, eⱼ₂, aⱼ₂, Iⱼ₂)
        Eⱼ₂ = AstroCoords._step3(fⱼ₂, eⱼ₂)
        rⱼ₂, νⱼ₂, Lⱼ₂ = AstroCoords._step4(aⱼ₂, eⱼ₂, η, Eⱼ₂, Lⱼ₂)
        _, δh = AstroCoords._step5(aⱼ₂, γ, γ′, θ, rⱼ₂, η, eⱼ₂, gⱼ₂, νⱼ₂, Lⱼ₂)

        f_fd, df_fd = value_and_jacobian(
            (x) ->
                collect(AstroCoords._step6(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8])),
            AutoFiniteDiff(),
            [γ′, θ, νⱼ₂, Lⱼ₂, eⱼ₂, gⱼ₂, hⱼ₂, δh],
        )

        f_ad, df_ad = value_and_jacobian(
            (x) ->
                collect(AstroCoords._step6(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8])),
            AutoMooncake(; config=nothing),
            [γ′, θ, νⱼ₂, Lⱼ₂, eⱼ₂, gⱼ₂, hⱼ₂, δh],
        )

        @test f_fd == f_ad
        @test df_fd ≈ df_ad atol = 1e-2
    end

    @testset "Step 7: Eccentricity Corrections" begin
        aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, fⱼ₂, Lⱼ₂ = AstroCoords._step1(state, μ)
        _, η, γ, γ′, θ = AstroCoords._step2(J2, Req, eⱼ₂, aⱼ₂, Iⱼ₂)
        Eⱼ₂ = AstroCoords._step3(fⱼ₂, eⱼ₂)
        rⱼ₂, νⱼ₂, Lⱼ₂ = AstroCoords._step4(aⱼ₂, eⱼ₂, η, Eⱼ₂, Lⱼ₂)

        f_fd, df_fd = value_and_jacobian(
            (x) -> collect(
                AstroCoords._step7(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9])
            ),
            AutoFiniteDiff(),
            [νⱼ₂, eⱼ₂, η, aⱼ₂, rⱼ₂, γ, γ′, θ, gⱼ₂],
        )

        f_ad, df_ad = value_and_jacobian(
            (x) -> collect(
                AstroCoords._step7(x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9])
            ),
            AutoMooncake(; config=nothing),
            [νⱼ₂, eⱼ₂, η, aⱼ₂, rⱼ₂, γ, γ′, θ, gⱼ₂],
        )

        @test f_fd == f_ad
        @test df_fd ≈ df_ad atol = 1e-2
    end

    @testset "Step 8: Inclination Corrections" begin
        aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, fⱼ₂, Lⱼ₂ = AstroCoords._step1(state, μ)
        _, η, γ, γ′, θ = AstroCoords._step2(J2, Req, eⱼ₂, aⱼ₂, Iⱼ₂)
        Eⱼ₂ = AstroCoords._step3(fⱼ₂, eⱼ₂)
        rⱼ₂, νⱼ₂, Lⱼ₂ = AstroCoords._step4(aⱼ₂, eⱼ₂, η, Eⱼ₂, Lⱼ₂)
        _, δh = AstroCoords._step5(aⱼ₂, γ, γ′, θ, rⱼ₂, η, eⱼ₂, gⱼ₂, νⱼ₂, Lⱼ₂)

        f_fd, df_fd = value_and_jacobian(
            (x) -> collect(AstroCoords._step8(x[1], x[2], x[3], x[4], x[5], x[6], x[7])),
            AutoFiniteDiff(),
            [γ′, θ, νⱼ₂, eⱼ₂, gⱼ₂, δh, Iⱼ₂],
        )

        f_ad, df_ad = value_and_jacobian(
            (x) -> collect(AstroCoords._step8(x[1], x[2], x[3], x[4], x[5], x[6], x[7])),
            AutoMooncake(; config=nothing),
            [γ′, θ, νⱼ₂, eⱼ₂, gⱼ₂, δh, Iⱼ₂],
        )

        @test f_fd == f_ad
        @test df_fd ≈ df_ad atol = 1e-2
    end

    @testset "Step 9: Final IOE Elements" begin
        aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, fⱼ₂, Lⱼ₂ = AstroCoords._step1(state, μ)
        _, η, γ, γ′, θ = AstroCoords._step2(J2, Req, eⱼ₂, aⱼ₂, Iⱼ₂)
        Eⱼ₂ = AstroCoords._step3(fⱼ₂, eⱼ₂)
        rⱼ₂, νⱼ₂, Lⱼ₂ = AstroCoords._step4(aⱼ₂, eⱼ₂, η, Eⱼ₂, Lⱼ₂)
        a, δh = AstroCoords._step5(aⱼ₂, γ, γ′, θ, rⱼ₂, η, eⱼ₂, gⱼ₂, νⱼ₂, Lⱼ₂)
        Σlgh = AstroCoords._step6(γ′, θ, νⱼ₂, Lⱼ₂, eⱼ₂, gⱼ₂, hⱼ₂, δh)
        _, _, _, _, δe, e″δL = AstroCoords._step7(νⱼ₂, eⱼ₂, η, aⱼ₂, rⱼ₂, γ, γ′, θ, gⱼ₂)
        δI, sin_half_I″_δh = AstroCoords._step8(γ′, θ, νⱼ₂, eⱼ₂, gⱼ₂, δh, Iⱼ₂)

        f_fd, df_fd = value_and_jacobian(
            (x) -> collect(
                AstroCoords._step9(
                    x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10]
                ),
            ),
            AutoFiniteDiff(),
            [a, eⱼ₂, δe, e″δL, Lⱼ₂, Iⱼ₂, δI, sin_half_I″_δh, hⱼ₂, Σlgh],
        )

        f_ad, df_ad = value_and_jacobian(
            (x) -> collect(
                AstroCoords._step9(
                    x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10]
                ),
            ),
            AutoMooncake(; config=nothing),
            [a, eⱼ₂, δe, e″δL, Lⱼ₂, Iⱼ₂, δI, sin_half_I″_δh, hⱼ₂, Σlgh],
        )

        @test f_fd == f_ad
        @test df_fd ≈ df_ad atol = 1e-2
    end
end

@testset "EDromo Transformation Differentiation" begin
    state = [
        -1076.225324679696,
        -6765.896364327722,
        -332.3087833503755,
        9.356857417032581,
        -3.3123476319597557,
        -1.1880157328553503,
    ]
    μ = 3.986004415e5
    cart_state = Cartesian(state)
    p = [0.0, 0.0]

    @testset "Differentiate wrt State" begin
        for backend in _BACKENDS
            @testset "$(backend[1]) Backend" begin
                for flag_time in (PhysicalTime(), ConstantTime(), LinearTime())
                    # Precompute for reverse pass
                    edromo_params = set_edromo_configurations(state, μ; flag_time=flag_time)
                    edromo_state_vec = Array(
                        params(EDromo(cart_state, μ; edromo_params...))
                    )

                    # Forward pass (Cartesian -> EDromo), including parameter calculation
                    to_edromo(x) = Array(
                        params(
                            EDromo(
                                Cartesian(x),
                                μ;
                                set_edromo_configurations(x, μ; flag_time=flag_time)...,
                            ),
                        ),
                    )

                    f_ad, df_ad = value_and_jacobian(x -> to_edromo(x), backend[2], state)
                    f_fd, df_fd = value_and_jacobian(
                        x -> to_edromo(x), AutoFiniteDiff(), state
                    )

                    @test f_ad ≈ f_fd rtol = 1e-8
                    @test df_ad ≈ df_fd rtol = 1e-6

                    # Reverse pass (EDromo -> Cartesian)
                    from_edromo(x) = Array(
                        params(Cartesian(EDromo(x...), μ; edromo_params...))
                    )

                    f_ad_rev, df_ad_rev = value_and_jacobian(
                        x -> from_edromo(x), backend[2], edromo_state_vec
                    )
                    f_fd_rev, df_fd_rev = value_and_jacobian(
                        x -> from_edromo(x), AutoFiniteDiff(), edromo_state_vec
                    )

                    @test f_ad_rev ≈ f_fd_rev rtol = 1e-8
                    @test df_ad_rev ≈ df_fd_rev rtol = 1e-6
                end
            end
        end
    end

    @testset "Differentiate wrt State Zygote" begin
        for flag_time in (PhysicalTime(), ConstantTime(), LinearTime())
            # Precompute for reverse pass
            edromo_params = set_edromo_configurations(state, μ; flag_time=flag_time)
            edromo_state_vec = Array(params(EDromo(cart_state, μ; edromo_params...)))

            # Forward pass (Cartesian -> EDromo), including parameter calculation
            to_edromo(x) = Array(
                params(
                    EDromo(
                        Cartesian(x),
                        μ;
                        set_edromo_configurations(x, μ; flag_time=flag_time)...,
                    ),
                ),
            )

            f_ad, df_ad = value_and_jacobian(x -> to_edromo(x), AutoZygote(), state)
            f_fd, df_fd = value_and_jacobian(x -> to_edromo(x), AutoFiniteDiff(), state)

            @test f_ad ≈ f_fd rtol = 1e-8
            @test df_ad ≈ df_fd rtol = 1e-6

            # Reverse pass (EDromo -> Cartesian)
            from_edromo(x) = Array(params(Cartesian(EDromo(x), μ; edromo_params...)))

            f_ad_rev, df_ad_rev = value_and_jacobian(
                x -> from_edromo(x), AutoZygote(), edromo_state_vec
            )
            f_fd_rev, df_fd_rev = value_and_jacobian(
                x -> from_edromo(x), AutoFiniteDiff(), edromo_state_vec
            )

            @test f_ad_rev ≈ f_fd_rev rtol = 1e-8
            @test df_ad_rev ≈ df_fd_rev rtol = 1e-6
        end
    end

    @testset "Differentiate wrt EDromo Parameters" begin
        for backend in _BACKENDS
            if backend[1] == "Enzyme"
                backend = ("Enzyme", AutoEnzyme(; function_annotation=Enzyme.Duplicated))
            end
            @testset "$(backend[1]) Backend" begin
                for flag_time in (PhysicalTime(), ConstantTime(), LinearTime())
                    # Differentiate wrt W and t₀
                    f_ad, df_ad = value_and_jacobian(
                        p -> Array(
                            params(
                                EDromo(
                                    cart_state,
                                    μ;
                                    set_edromo_configurations(
                                        state, μ; W=p[1], t₀=p[2], flag_time=flag_time
                                    )...,
                                ),
                            ),
                        ),
                        backend[2],
                        p,
                    )
                    f_fd, df_fd = value_and_jacobian(
                        p -> Array(
                            params(
                                EDromo(
                                    cart_state,
                                    μ;
                                    set_edromo_configurations(
                                        state, μ; W=p[1], t₀=p[2], flag_time=flag_time
                                    )...,
                                ),
                            ),
                        ),
                        AutoFiniteDiff(),
                        p,
                    )

                    @test f_ad ≈ f_fd rtol = 1e-8
                    @test df_ad ≈ df_fd rtol = 1e-6
                end
            end
        end
    end

    @testset "Differentiate wrt EDromo Parameters Zygote" begin
        for flag_time in (PhysicalTime(), ConstantTime(), LinearTime())
            # Differentiate wrt W and t₀
            to_edromo_params(p) = params(
                EDromo(
                    cart_state,
                    μ;
                    set_edromo_configurations(
                        state, μ; W=p[1], t₀=p[2], flag_time=flag_time
                    )...,
                ),
            )

            f_ad, df_ad = value_and_jacobian(p -> to_edromo_params(p), AutoZygote(), p)
            f_fd, df_fd = value_and_jacobian(p -> to_edromo_params(p), AutoFiniteDiff(), p)

            @test f_ad ≈ f_fd rtol = 1e-8
            @test df_ad ≈ df_fd rtol = 1e-6
        end
    end

    @testset "Differentiate wrt μ" begin
        for backend in _BACKENDS
            if backend[1] == "Enzyme"
                backend = ("Enzyme", AutoEnzyme(; function_annotation=Enzyme.Duplicated))
            end
            @testset "$(backend[1]) Backend" begin
                for flag_time in (PhysicalTime(), ConstantTime(), LinearTime())
                    # Setup for reverse pass
                    edromo_params = set_edromo_configurations(state, μ; flag_time=flag_time)
                    edromo_state = EDromo(cart_state, μ; edromo_params...)

                    f_ad, df_ad = value_and_derivative(
                        m -> Array(
                            params(
                                EDromo(
                                    cart_state,
                                    m;
                                    set_edromo_configurations(
                                        state, m; flag_time=flag_time
                                    )...,
                                ),
                            ),
                        ),
                        backend[2],
                        μ,
                    )
                    f_fd, df_fd = value_and_derivative(
                        m -> Array(
                            params(
                                EDromo(
                                    cart_state,
                                    m;
                                    set_edromo_configurations(
                                        state, m; flag_time=flag_time
                                    )...,
                                ),
                            ),
                        ),
                        AutoFiniteDiff(),
                        μ,
                    )

                    @test f_ad ≈ f_fd rtol = 1e-8
                    @test df_ad ≈ df_fd rtol = 1e-3

                    f_ad_rev, df_ad_rev = value_and_derivative(
                        m -> Array(params(Cartesian(edromo_state, m; edromo_params...))),
                        backend[2],
                        μ,
                    )
                    f_fd_rev, df_fd_rev = value_and_derivative(
                        m -> Array(params(Cartesian(edromo_state, m; edromo_params...))),
                        AutoFiniteDiff(),
                        μ,
                    )

                    @test f_ad_rev ≈ f_fd_rev rtol = 1e-8
                    @test df_ad_rev ≈ df_fd_rev rtol = 1e-3
                end
            end
        end
    end

    @testset "Differentiate wrt μ Zygote" begin
        for flag_time in (PhysicalTime(), ConstantTime(), LinearTime())
            # Setup for reverse pass
            edromo_params = set_edromo_configurations(state, μ; flag_time=flag_time)
            edromo_state = EDromo(cart_state, μ; edromo_params...)

            # Forward pass (Cartesian -> EDromo)
            to_edromo_μ(m) = params(
                EDromo(
                    cart_state,
                    m;
                    set_edromo_configurations(state, m; flag_time=flag_time)...,
                ),
            )

            f_ad, df_ad = value_and_derivative(m -> to_edromo_μ(m), AutoZygote(), μ)
            f_fd, df_fd = value_and_derivative(m -> to_edromo_μ(m), AutoFiniteDiff(), μ)

            @test f_ad ≈ f_fd rtol = 1e-8
            @test df_ad ≈ df_fd atol = 1e-3

            # Reverse pass (EDromo -> Cartesian)
            from_edromo_μ(m) = params(Cartesian(edromo_state, m; edromo_params...))

            try
                f_ad_rev, df_ad_rev = value_and_derivative(
                    m -> from_edromo_μ(m), AutoZygote(), μ
                )
                f_fd_rev, df_fd_rev = value_and_derivative(
                    m -> from_edromo_μ(m), AutoFiniteDiff(), μ
                )

                @test f_ad_rev ≈ f_fd_rev rtol = 1e-8
                @test df_ad_rev ≈ df_fd_rev atol = 1e-3
            catch err
                @test err isa MethodError
                @test startswith(
                    sprint(showerror, err),
                    "MethodError: no method matching iterate(::Nothing)",
                )
            end
        end
    end
end
