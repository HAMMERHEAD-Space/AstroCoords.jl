@testset "Keplerian Coordinate Set" begin
    @testset "Inner constructor Keplerian{T}" begin
        k = Keplerian{Float64}(7000.0, 0.01, 0.1, 0.2, 0.3, 0.4)
        @test k isa Keplerian{Float64}
        @test k.a ≈ 7000.0
        @test k.e ≈ 0.01
        @test k.i ≈ 0.1
        @test k.Ω ≈ 0.2
        @test k.ω ≈ 0.3
        @test k.f ≈ 0.4
    end

    @testset "Constructor with type promotion" begin
        # Float32 + Float64 should promote to Float64
        k = Keplerian(7000.0f0, 0.01, 0.1f0, 0.2, 0.5f0, 0.4)
        @test k isa Keplerian{Float64}
        @test k.a ≈ 7000.0
        @test k.e ≈ 0.01
        @test k.i ≈ 0.1
        @test k.Ω ≈ 0.2
        @test k.ω ≈ 0.5
        @test k.f ≈ 0.4
    end

    @testset "StaticVector constructor" begin
        vec = SVector{6}(7000.0, 0.01, 0.1, 0.2, 0.3, 0.4)
        k = Keplerian(vec)
        @test k isa Keplerian{Float64}
        @test k.a ≈ 7000.0
        @test k.e ≈ 0.01
        @test k.i ≈ 0.1
        @test k.Ω ≈ 0.2
        @test k.ω ≈ 0.3
        @test k.f ≈ 0.4
    end

    @testset "Constructor from vector" begin
        vec = [7000.0, 0.01, 0.1, 0.2, 0.3, 0.4]
        k = Keplerian(vec)
        @test k isa Keplerian{Float64}
        @test k.a ≈ 7000.0
        @test k.e ≈ 0.01
        @test k.i ≈ 0.1
        @test k.Ω ≈ 0.2
        @test k.ω ≈ 0.3
        @test k.f ≈ 0.4
    end

    @testset "Base.one(Keplerian)" begin
        k_one = Base.one(Keplerian)
        @test k_one isa Keplerian{Float64}
        @test k_one.a ≈ 0.0
        @test k_one.e ≈ 0.0
        @test k_one.i ≈ 0.0
        @test k_one.Ω ≈ 0.0
        @test k_one.ω ≈ 0.0
        @test k_one.f ≈ 0.0

        k_one_f32 = Base.one(Keplerian; T=Float32)
        @test k_one_f32 isa Keplerian{Float32}
        @test k_one_f32.a ≈ 0.0f0
    end

    @testset "Base.getindex for all 6 indices" begin
        k = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.4)

        @test k[1] ≈ 7000.0  # a
        @test k[2] ≈ 0.01    # e
        @test k[3] ≈ 0.1     # i
        @test k[4] ≈ 0.2     # Ω
        @test k[5] ≈ 0.3     # ω
        @test k[6] ≈ 0.4     # f

        @test_throws BoundsError k[0]
        @test_throws BoundsError k[7]
    end

    @testset "Type stability" begin
        k_f64 = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.4)
        @test typeof(k_f64.a) === Float64
        @test typeof(k_f64.e) === Float64

        k_f32 = Keplerian(7000.0f0, 0.01f0, 0.1f0, 0.2f0, 0.3f0, 0.4f0)
        @test typeof(k_f32.a) === Float32
        @test typeof(k_f32.e) === Float32
    end
end
