@testset "Round Trip Coordinate Changes" begin
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

    for coord in filter(T -> T ∉ (EDromo, KustaanheimoStiefel), AstroCoords.COORD_TYPES)
        coord_state = coord(cart_state, μ)
        for coord2 in
            filter(T -> T ∉ (EDromo, KustaanheimoStiefel), AstroCoords.COORD_TYPES)
            coord_state2 = coord2(coord_state, μ)
            coord_state_round_trip2 = coord(coord_state2, μ)
            @test params(coord_state) ≈ params(coord_state_round_trip2) rtol = 1e-14
        end
    end
end
