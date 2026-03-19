for set in _COORDINATE_SETS
    @test length(check_allocs(set, (Vector{Number},))) == 0
end

@testset "Transformation Allocations" begin
    # Test over all coordinate sets
    for set1 in filter(
        T -> T ∉ (EDromo, KustaanheimoStiefel, StiefelScheifele, GEqOE),
        AstroCoords.COORD_TYPES,
    )
        for set2 in filter(
            T -> T ∉ (EDromo, KustaanheimoStiefel, StiefelScheifele, GEqOE),
            AstroCoords.COORD_TYPES,
        )
            @test length(check_allocs(set1, (set2{Float64}, Float64))) == 0
        end
    end
end

@testset "Quantity Allocations" begin
    for set in filter(
        T -> T ∉ (EDromo, KustaanheimoStiefel, StiefelScheifele, GEqOE),
        AstroCoords.COORD_TYPES,
    )
        @test length(check_allocs(meanMotion, (set{Float64}, Float64))) == 0
        @test length(check_allocs(orbitalPeriod, (set{Float64}, Float64))) == 0
        @test length(check_allocs(orbitalNRG, (set{Float64}, Float64))) == 0
        @test length(check_allocs(angularMomentumVector, (set{Float64}, Float64))) == 0
        @test length(check_allocs(angularMomentumQuantity, (set{Float64}, Float64))) == 0
    end
end

@testset "EDromo Allocs" begin
    state = Cartesian([
        -1076.225324679696,
        -6765.896364327722,
        -332.3087833503755,
        9.356857417032581,
        -3.3123476319597557,
        -1.1880157328553503,
    ])
    μ = 3.986004415e5
    edromo_params = RegularizedCoordinateConfig(params(state), μ)
    # Compute phi separately
    ϕ = compute_initial_phi(params(state), μ, edromo_params)
    edromo_state = EDromo(state, μ, ϕ, edromo_params)

    @test length(
        check_allocs(EDromo, (Cartesian{Float64}, Float64, Float64, typeof(edromo_params)))
    ) == 0
    @test length(
        check_allocs(Cartesian, (EDromo{Float64}, Float64, Float64, typeof(edromo_params)))
    ) == 0
    @test length(
        check_allocs(meanMotion, (EDromo{Float64}, Float64, Float64, typeof(edromo_params)))
    ) == 0
    @test length(
        check_allocs(
            orbitalPeriod, (EDromo{Float64}, Float64, Float64, typeof(edromo_params))
        ),
    ) == 0
    @test length(
        check_allocs(orbitalNRG, (EDromo{Float64}, Float64, Float64, typeof(edromo_params)))
    ) == 0
    @test length(
        check_allocs(
            angularMomentumVector,
            (EDromo{Float64}, Float64, Float64, typeof(edromo_params)),
        ),
    ) == 0
    @test length(
        check_allocs(
            angularMomentumQuantity,
            (EDromo{Float64}, Float64, Float64, typeof(edromo_params)),
        ),
    ) == 0

    # Test allocation for config construction
    config_constructor(state, μ) = RegularizedCoordinateConfig(state, μ)
    @test length(check_allocs(config_constructor, (typeof(state), typeof(μ)))) == 0
end

@testset "Kustaanheimo-Stiefel Allocs" begin
    state = Cartesian([
        -1076.225324679696,
        -6765.896364327722,
        -332.3087833503755,
        9.356857417032581,
        -3.3123476319597557,
        -1.1880157328553503,
    ])
    μ = 3.986004415e5
    ks_params = RegularizedCoordinateConfig(state, μ)
    ks_state = KustaanheimoStiefel(state, μ, ks_params)

    @test length(
        check_allocs(KustaanheimoStiefel, (Cartesian{Float64}, Float64, typeof(ks_params)))
    ) == 0
    @test length(
        check_allocs(Cartesian, (KustaanheimoStiefel{Float64}, Float64, typeof(ks_params)))
    ) == 0
    @test length(
        check_allocs(meanMotion, (KustaanheimoStiefel{Float64}, Float64, typeof(ks_params)))
    ) == 0
    @test length(
        check_allocs(
            orbitalPeriod, (KustaanheimoStiefel{Float64}, Float64, typeof(ks_params))
        ),
    ) == 0
    @test length(
        check_allocs(orbitalNRG, (KustaanheimoStiefel{Float64}, Float64, typeof(ks_params)))
    ) == 0
    @test length(
        check_allocs(
            angularMomentumVector,
            (KustaanheimoStiefel{Float64}, Float64, typeof(ks_params)),
        ),
    ) == 0
    @test length(
        check_allocs(
            angularMomentumQuantity,
            (KustaanheimoStiefel{Float64}, Float64, typeof(ks_params)),
        ),
    ) == 0

    # Test allocation for config construction 
    config_constructor_ks(state, μ) = RegularizedCoordinateConfig(state, μ)
    @test length(check_allocs(config_constructor_ks, (typeof(state), typeof(μ)))) == 0
end

@testset "Stiefel-Scheifele Allocs" begin
    state = Cartesian([
        -1076.225324679696,
        -6765.896364327722,
        -332.3087833503755,
        9.356857417032581,
        -3.3123476319597557,
        -1.1880157328553503,
    ])
    μ = 3.986004415e5
    ss_params = RegularizedCoordinateConfig(params(state), μ)
    # Compute phi separately
    ϕ = compute_initial_phi(params(state), μ, ss_params)
    ss_state = StiefelScheifele(state, μ, ϕ, ss_params)

    @test length(
        check_allocs(
            StiefelScheifele, (Cartesian{Float64}, Float64, Float64, typeof(ss_params))
        ),
    ) == 0
    @test length(
        check_allocs(
            Cartesian, (StiefelScheifele{Float64}, Float64, Float64, typeof(ss_params))
        ),
    ) == 0
    @test length(
        check_allocs(
            meanMotion, (StiefelScheifele{Float64}, Float64, Float64, typeof(ss_params))
        ),
    ) == 0
    @test length(
        check_allocs(
            orbitalPeriod, (StiefelScheifele{Float64}, Float64, Float64, typeof(ss_params))
        ),
    ) == 0
    @test length(
        check_allocs(
            orbitalNRG, (StiefelScheifele{Float64}, Float64, Float64, typeof(ss_params))
        ),
    ) == 0
    @test length(
        check_allocs(
            angularMomentumVector,
            (StiefelScheifele{Float64}, Float64, Float64, typeof(ss_params)),
        ),
    ) == 0
    @test length(
        check_allocs(
            angularMomentumQuantity,
            (StiefelScheifele{Float64}, Float64, Float64, typeof(ss_params)),
        ),
    ) == 0

    # Test allocation for config construction
    config_constructor_ss(state, μ) = RegularizedCoordinateConfig(state, μ)
    @test length(check_allocs(config_constructor_ss, (typeof(state), typeof(μ)))) == 0
end

@testset "OrbitState Allocs" begin
    u0 = Cartesian([
        -1076.225324679696,
        -6765.896364327722,
        -332.3087833503755,
        9.356857417032581,
        -3.3123476319597557,
        -1.1880157328553503,
    ])
    μ = 3.986004415e5
    ep = 0.0
    s = OrbitState(u0, ep, :ICRF)

    # Construction: Symbol convenience path
    @test length(
        check_allocs((c, e) -> OrbitState(c, e, :ICRF), (typeof(u0), typeof(ep)))
    ) == 0

    # convert_coords: Cartesian → Keplerian
    @test length(check_allocs(convert_coords, (typeof(s), Type{Keplerian}, typeof(μ)))) == 0

    # convert_coords: Keplerian → Cartesian
    s_kep = convert_coords(s, Keplerian, μ)
    @test length(
        check_allocs(convert_coords, (typeof(s_kep), Type{Cartesian}, typeof(μ)))
    ) == 0

    # Accessors
    @test length(check_allocs(coords, (typeof(s),))) == 0
    @test length(check_allocs(epoch, (typeof(s),))) == 0
    @test length(check_allocs(frame, (typeof(s),))) == 0

    frames = FrameSystem{2,Float64}()
    add_axes!(frames, :ICRF, 1)
    add_axes_fixed_angles!(frames, :ITRF, 2, 1, [0.0, 0.0, 0.0], :ZYX)
    cr_fwd = compile_rotation6(frames, :ICRF, :ITRF)
    cr_inv = compile_rotation6(frames, :ITRF, :ICRF)

    @test length(
        check_allocs(change_frame, (typeof(s), Val{:ITRF}, typeof(cr_fwd), typeof(μ)))
    ) == 0

    s_itrf = change_frame(s, Val{:ITRF}(), cr_fwd, μ)
    @test length(
        check_allocs(change_frame, (typeof(s_itrf), Val{:ICRF}, typeof(cr_inv), typeof(μ)))
    ) == 0
end

@testset "GEqOE Allocs" begin
    state = Cartesian([
        -1076.225324679696,
        -6765.896364327722,
        -332.3087833503755,
        9.356857417032581,
        -3.3123476319597557,
        -1.1880157328553503,
    ])
    μ = 3.986004415e5
    geqoe_config = RegularizedCoordinateConfig(; W=0.0)

    geqoe_state = GEqOE(state, μ, geqoe_config)

    @test length(
        check_allocs(GEqOE, (Cartesian{Float64}, Float64, typeof(geqoe_config)))
    ) == 0
    @test length(
        check_allocs(Cartesian, (GEqOE{Float64}, Float64, typeof(geqoe_config)))
    ) == 0
    @test length(
        check_allocs(meanMotion, (GEqOE{Float64}, Float64, typeof(geqoe_config)))
    ) == 0
    @test length(
        check_allocs(orbitalPeriod, (GEqOE{Float64}, Float64, typeof(geqoe_config)))
    ) == 0
    @test length(
        check_allocs(orbitalNRG, (GEqOE{Float64}, Float64, typeof(geqoe_config)))
    ) == 0
    @test length(
        check_allocs(angularMomentumVector, (GEqOE{Float64}, Float64, typeof(geqoe_config)))
    ) == 0
    @test length(
        check_allocs(
            angularMomentumQuantity, (GEqOE{Float64}, Float64, typeof(geqoe_config))
        ),
    ) == 0
end
