@testset "Stiefel-Scheifele" begin
    @testset "Round Trip" begin
        # Initial state
        state = Cartesian([
            -1076.225324679696,
            -6765.896364327722,
            -332.3087833503755,
            9.356857417032581,
            -3.3123476319597557,
            -1.1880157328553503,
        ])
        μ = 3.986004415e5

        # Get parameters
        ss_params = set_stiefelscheifele_configurations(state, μ)

        # Forward transformation
        ss_state = StiefelScheifele(state, μ; ss_params...)

        # Backward transformation
        final_state = Cartesian(ss_state, μ; ss_params...)

        @test final_state ≈ state
    end

    @testset "Configurations" begin
        state = Cartesian([1000.0, 5000.0, -3000.0, -5.0, 2.0, 7.0])
        μ = 3.986004415e5
        r_norm = norm(state[1:3])
        
        configs = set_stiefelscheifele_configurations(state, μ)
        
        @test configs.DU == r_norm
        @test configs.TU == 1 / sqrt(μ / r_norm^3)
        @test configs.W == 0.0
        @test configs.ϕ == 0.0
        @test configs.t₀ == 0.0
        @test configs.flag_time == 0
    end
end 