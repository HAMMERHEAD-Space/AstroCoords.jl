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
    @test NRG ≈ -8.146487030627135 rtol=1e-14
    @test h_vec ≈ [6937.269116080104; -4387.938522053871; 66872.35998509153] rtol=1e-14
    @test h ≈ 67374.26984567562 rtol=1e-14

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

        @test NRG ≈ NRG2 rtol=1e-14
        @test h_vec ≈ h_vec2 rtol=1e-14
        @test h ≈ h rtol=1e-14
    end
end

