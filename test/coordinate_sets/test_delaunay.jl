@testset "Delaunay Coordinate Set" begin
    @testset "Inner constructor Delaunay{T}" begin
        d = Delaunay{Float64}(1.0, 2.0, 3.0, 0.5, 1.0, 2.0)
        @test d isa Delaunay{Float64}
        @test d.L ≈ 1.0
        @test d.G ≈ 2.0
        @test d.H ≈ 3.0
        @test d.M ≈ 0.5
        @test d.ω ≈ 1.0
        @test d.Ω ≈ 2.0
    end

    @testset "Constructor with type promotion" begin
        # Float32 + Float64 should promote to Float64
        d = Delaunay(1.0f0, 2.0, 3.0f0, 0.5, 1.0f0, 2.0)
        @test d isa Delaunay{Float64}
        @test d.L ≈ 1.0
        @test d.G ≈ 2.0
        @test d.H ≈ 3.0
        @test d.M ≈ 0.5
        @test d.ω ≈ 1.0
        @test d.Ω ≈ 2.0
    end

    @testset "StaticVector constructor" begin
        vec = SVector{6}(1.0, 2.0, 3.0, 0.5, 1.0, 2.0)
        d = Delaunay(vec)
        @test d isa Delaunay{Float64}
        @test d.L ≈ 1.0
        @test d.G ≈ 2.0
        @test d.H ≈ 3.0
        @test d.M ≈ 0.5
        @test d.ω ≈ 1.0
        @test d.Ω ≈ 2.0
    end

    @testset "Constructor from vector" begin
        vec = [1.0, 2.0, 3.0, 0.5, 1.0, 2.0]
        d = Delaunay(vec)
        @test d isa Delaunay{Float64}
        @test d.L ≈ 1.0
        @test d.G ≈ 2.0
        @test d.H ≈ 3.0
        @test d.M ≈ 0.5
        @test d.ω ≈ 1.0
        @test d.Ω ≈ 2.0
    end

    @testset "Base.one(Delaunay)" begin
        d_one = Base.one(Delaunay)
        @test d_one isa Delaunay{Float64}
        @test d_one.L ≈ 0.0
        @test d_one.G ≈ 0.0
        @test d_one.H ≈ 0.0
        @test d_one.M ≈ 0.0
        @test d_one.ω ≈ 0.0
        @test d_one.Ω ≈ 0.0

        d_one_f32 = Base.one(Delaunay; T=Float32)
        @test d_one_f32 isa Delaunay{Float32}
        @test d_one_f32.L ≈ 0.0f0
    end

    @testset "Base.getindex for all 6 indices" begin
        d = Delaunay(1.0, 2.0, 3.0, 0.5, 1.0, 2.0)

        @test d[1] ≈ 1.0  # L
        @test d[2] ≈ 2.0  # G
        @test d[3] ≈ 3.0  # H
        @test d[4] ≈ 0.5  # M
        @test d[5] ≈ 1.0  # ω
        @test d[6] ≈ 2.0  # Ω

        @test_throws BoundsError d[0]
        @test_throws BoundsError d[7]
    end
end
