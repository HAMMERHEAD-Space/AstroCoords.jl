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

    for coord in _COORDINATE_SETS
        println(coord)
        coord_state = coord(cart_state, μ)
        cart_state_round_trip = Cartesian(coord_state, μ)
        @test params(cart_state) ≈ params(cart_state_round_trip) rtol=1e-14
    end
end
