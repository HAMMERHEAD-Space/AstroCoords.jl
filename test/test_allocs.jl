for set in _COORDINATE_SETS
    @test length(check_allocs(set, (Vector{Number},))) == 0
end

@testset "Transformation Allocations" begin
    # Test over all coordinate sets
    for set1 in filter(
        T -> T ∉ (EDromo, KustaanheimoStiefel, StiefelScheifele), AstroCoords.COORD_TYPES
    )
        for set2 in filter(
            T -> T ∉ (EDromo, KustaanheimoStiefel, StiefelScheifele),
            AstroCoords.COORD_TYPES,
        )
            @test length(check_allocs(set1, (set2{Float64}, Float64))) == 0
        end
    end
end

@testset "Quantity Allocations" begin
    for set in filter(
        T -> T ∉ (EDromo, KustaanheimoStiefel, StiefelScheifele), AstroCoords.COORD_TYPES
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
