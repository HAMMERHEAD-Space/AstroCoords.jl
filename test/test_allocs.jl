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
    edromo_params = set_edromo_configurations(state, μ)
    edromo_state = EDromo(state, μ; edromo_params...)

    to_edromo(x, μ) = EDromo(x, μ; edromo_params...)
    from_edromo(x, μ) = Cartesian(x, μ; edromo_params...)
    mm_edromo(x, μ) = meanMotion(x, μ; edromo_params...)
    op_edromo(x, μ) = orbitalPeriod(x, μ; edromo_params...)
    on_edromo(x, μ) = orbitalNRG(x, μ; edromo_params...)
    av_edromo(x, μ) = angularMomentumVector(x, μ; edromo_params...)
    aq_edromo(x, μ) = angularMomentumQuantity(x, μ; edromo_params...)

    @test length(check_allocs(to_edromo, (Cartesian{Float64}, Float64))) == 0
    @test length(check_allocs(from_edromo, (EDromo{Float64}, Float64))) == 0
    @test length(check_allocs(mm_edromo, (EDromo{Float64}, Float64))) == 0
    @test length(check_allocs(op_edromo, (EDromo{Float64}, Float64))) == 0
    @test length(check_allocs(on_edromo, (EDromo{Float64}, Float64))) == 0
    @test length(check_allocs(av_edromo, (EDromo{Float64}, Float64))) == 0
    @test length(check_allocs(aq_edromo, (EDromo{Float64}, Float64))) == 0
    @test length(check_allocs(set_edromo_configurations, (typeof(state), typeof(μ)))) == 0
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
    ks_params = set_ks_configurations(state, μ)
    ks_state = KustaanheimoStiefel(state, μ; ks_params...)

    to_ks(x, μ) = KustaanheimoStiefel(x, μ; ks_params...)
    from_ks(x, μ) = Cartesian(x, μ; ks_params...)
    mm_ks(x, μ) = meanMotion(x, μ; ks_params...)
    op_ks(x, μ) = orbitalPeriod(x, μ; ks_params...)
    on_ks(x, μ) = orbitalNRG(x, μ; ks_params...)
    av_ks(x, μ) = angularMomentumVector(x, μ; ks_params...)
    aq_ks(x, μ) = angularMomentumQuantity(x, μ; ks_params...)

    @test length(check_allocs(to_ks, (Cartesian{Float64}, Float64))) == 0
    @test length(check_allocs(from_ks, (KustaanheimoStiefel{Float64}, Float64))) == 0
    @test length(check_allocs(mm_ks, (KustaanheimoStiefel{Float64}, Float64))) == 0
    @test length(check_allocs(op_ks, (KustaanheimoStiefel{Float64}, Float64))) == 0
    @test length(check_allocs(on_ks, (KustaanheimoStiefel{Float64}, Float64))) == 0
    @test length(check_allocs(av_ks, (KustaanheimoStiefel{Float64}, Float64))) == 0
    @test length(check_allocs(aq_ks, (KustaanheimoStiefel{Float64}, Float64))) == 0
    @test length(check_allocs(set_ks_configurations, (typeof(state), typeof(μ)))) == 0
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
    ss_params = set_stiefelscheifele_configurations(state, μ)
    ss_state = StiefelScheifele(state, μ; ss_params...)

    to_ss(x, μ) = StiefelScheifele(x, μ; ss_params...)
    from_ss(x, μ) = Cartesian(x, μ; ss_params...)
    mm_ss(x, μ) = meanMotion(x, μ; ss_params...)
    op_ss(x, μ) = orbitalPeriod(x, μ; ss_params...)
    on_ss(x, μ) = orbitalNRG(x, μ; ss_params...)
    av_ss(x, μ) = angularMomentumVector(x, μ; ss_params...)
    aq_ss(x, μ) = angularMomentumQuantity(x, μ; ss_params...)

    @test length(check_allocs(to_ss, (Cartesian{Float64}, Float64))) == 0
    @test length(check_allocs(from_ss, (StiefelScheifele{Float64}, Float64))) == 0
    @test length(check_allocs(mm_ss, (StiefelScheifele{Float64}, Float64))) == 0
    @test length(check_allocs(op_ss, (StiefelScheifele{Float64}, Float64))) == 0
    @test length(check_allocs(on_ss, (StiefelScheifele{Float64}, Float64))) == 0
    @test length(check_allocs(av_ss, (StiefelScheifele{Float64}, Float64))) == 0
    @test length(check_allocs(aq_ss, (StiefelScheifele{Float64}, Float64))) == 0
    @test length(
        check_allocs(set_stiefelscheifele_configurations, (typeof(state), typeof(μ)))
    ) == 0
end
