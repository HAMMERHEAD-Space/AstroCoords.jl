for set in _COORDINATE_SETS
    @test length(check_allocs(set, (Vector{Number},))) == 0
end

for set1 in filter(T -> T != EDromo, AstroCoords.COORD_TYPES)
    for set2 in filter(T -> T != EDromo, AstroCoords.COORD_TYPES)
        @test length(check_allocs(set1, (set2{Float64}, Float64))) == 0
    end
end

for set in filter(T -> T != EDromo, AstroCoords.COORD_TYPES)
    @test length(check_allocs(meanMotion, (set{Float64}, Float64))) == 0
    @test length(check_allocs(orbitalPeriod, (set{Float64}, Float64))) == 0
    @test length(check_allocs(orbitalNRG, (set{Float64}, Float64))) == 0
    @test length(check_allocs(angularMomentumVector, (set{Float64}, Float64))) == 0
    @test length(check_allocs(angularMomentumQuantity, (set{Float64}, Float64))) == 0
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
    edromo_params = get_edromo_defaults(state, μ)
    edromo_state = EDromo(state, μ; edromo_params...)

    to_edromo(x, μ) = EDromo(x, μ; edromo_params...)
    from_edromo(x, μ) = Cartesian(x, μ; edromo_params...)
    mm_edromo(x, μ) = meanMotion(x, μ; edromo_params...)
    op_edromo(x, μ) = orbitalPeriod(x, μ; edromo_params...)
    on_edromo(x, μ) = orbitalNRG(x, μ; edromo_params...)
    av_edromo(x, μ) = angularMomentumVector(x, μ; edromo_params...)
    aq_edromo(x, μ) = angularMomentumQuantity(x, μ; edromo_params...)

    @test length(check_allocs(get_edromo_defaults, (typeof(state), typeof(μ)))) == 0
    @test length(check_allocs(to_edromo, (typeof(state), typeof(μ)))) == 0
    @test length(check_allocs(from_edromo, (typeof(edromo_state), typeof(μ)))) == 0
    @test length(check_allocs(mm_edromo, (typeof(edromo_state), typeof(μ)))) == 0
    @test length(check_allocs(op_edromo, (typeof(edromo_state), typeof(μ)))) == 0
    @test length(check_allocs(on_edromo, (typeof(edromo_state), typeof(μ)))) == 0
    @test length(check_allocs(av_edromo, (typeof(edromo_state), typeof(μ)))) == 0
    @test length(check_allocs(aq_edromo, (typeof(edromo_state), typeof(μ)))) == 0
end
