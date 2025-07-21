@testset "Stiefel-Scheifele Round Trip" begin
    # Initial state
    base_state_vec = [
        -1076.225324679696,
        -6765.896364327722,
        -332.3087833503755,
        9.356857417032581,
        -3.3123476319597557,
        -1.1880157328553503,
    ]
    μ = 3.986004415e5
    base_cart_state = Cartesian(base_state_vec)

    # All coordinate types except EDromo itself
    coord_types_to_test = filter(
        T -> T ∉ (EDromo, KustaanheimoStiefel, StiefelScheifele), AstroCoords.COORD_TYPES
    )

    # Parameter sets to test
    time_flags = [PhysicalTime(), LinearTime()]
    W_values = [0.0, 1e-8]
    t0_values = [0.0, 100.0]

    for FromCoord in coord_types_to_test
        @testset "Roundtrip from $FromCoord" begin
            # Create the initial state for this coordinate type
            from_state = FromCoord(base_cart_state, μ)

            for flag in time_flags
                for W in W_values
                    for t₀ in t0_values
                        @testset "Params: time_flag=$(typeof(flag)), W=$W, t₀=$t₀" begin
                            # Get the full set of parameters for this test case
                            defaults = RegularizedCoordinateConfig(
                                base_state_vec, μ; W=W, t₀=t₀, flag_time=flag
                            )

                            # Perform the round trip
                            ss_state = StiefelScheifele(from_state, μ, defaults)
                            roundtrip_state = FromCoord(ss_state, μ, defaults)

                            # Test for numerical equality
                            @test params(roundtrip_state) ≈ params(from_state) atol=1e-8
                        end
                    end
                end
            end
        end
    end
end
