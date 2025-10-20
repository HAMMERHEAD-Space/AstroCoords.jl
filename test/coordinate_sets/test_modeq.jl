@testset "ModEq and J2EqOE Coordinate Sets" begin
    @testset "ModEq - Inner Constructor" begin
        modeq = ModEq{Float64}(7000.0, 0.1, 0.2, 0.05, 0.06, 1.0)
        @test modeq isa ModEq{Float64}
        @test modeq.p == 7000.0
        @test modeq.f == 0.1
        @test modeq.g == 0.2
        @test modeq.h == 0.05
        @test modeq.k == 0.06
        @test modeq.L == 1.0
    end

    @testset "ModEq - Constructor from Vector" begin
        vec = [7000.0, 0.1, 0.2, 0.05, 0.06, 1.0]
        modeq = ModEq(vec)
        @test modeq isa ModEq{Float64}
        @test modeq.p == 7000.0
        @test modeq.L == 1.0

        # Test with Int vector
        vec_int = [7000, 0, 0, 0, 0, 1]
        modeq_int = ModEq(vec_int)
        @test modeq_int isa ModEq{Int}
    end

    @testset "ModEq - Type Promotion" begin
        # Mixed types
        modeq = ModEq(7000.0, 0.1f0, 0.2, 0, 0.06, 1.0)
        @test modeq isa ModEq
        @test eltype(modeq) <: AbstractFloat

        # All Float32
        modeq_f32 = ModEq(7000.0f0, 0.1f0, 0.2f0, 0.05f0, 0.06f0, 1.0f0)
        @test modeq_f32 isa ModEq{Float32}
    end

    @testset "ModEq - StaticVector Constructor" begin
        svec = SVector{6}(7000.0, 0.1, 0.2, 0.05, 0.06, 1.0)
        modeq = ModEq(svec)
        @test modeq isa ModEq{Float64}
        @test modeq.p == 7000.0
    end

    @testset "ModEq - Constructor from Vector (line 38)" begin
        # Test the specific constructor that takes AbstractVector
        vec = [7000.0, 0.1, 0.2, 0.05, 0.06, 1.0]
        modeq = ModEq(vec)
        @test modeq isa ModEq{Float64}
        @test length(params(modeq)) == 6
    end

    @testset "ModEq - Base.one()" begin
        modeq_one = one(ModEq)
        @test modeq_one isa ModEq{Float64}
        @test all(params(modeq_one) .== 0.0)

        modeq_one_f32 = one(ModEq, T=Float32)
        @test modeq_one_f32 isa ModEq{Float32}
    end

    @testset "ModEq - params() Method" begin
        modeq = ModEq(7000.0, 0.1, 0.2, 0.05, 0.06, 1.0)
        p = params(modeq)

        @test p isa SVector{6,Float64}
        @test p[1] == 7000.0
        @test p[2] == 0.1
        @test p[3] == 0.2
        @test p[4] == 0.05
        @test p[5] == 0.06
        @test p[6] == 1.0
    end

    @testset "ModEq - Base.getindex" begin
        modeq = ModEq(7000.0, 0.1, 0.2, 0.05, 0.06, 1.0)

        @test modeq[1] == 7000.0  # p
        @test modeq[2] == 0.1     # f
        @test modeq[3] == 0.2     # g
        @test modeq[4] == 0.05    # h
        @test modeq[5] == 0.06    # k
        @test modeq[6] == 1.0     # L

        @test_throws BoundsError modeq[0]
        @test_throws BoundsError modeq[7]
    end

    @testset "J2EqOE - Inner Constructor" begin
        j2eq = J2EqOE{Float64}(0.001, 0.1, 0.2, 0.05, 0.06, 1.0)
        @test j2eq isa J2EqOE{Float64}
        @test j2eq.n == 0.001
        @test j2eq.h == 0.1
        @test j2eq.k == 0.2
        @test j2eq.p == 0.05
        @test j2eq.q == 0.06
        @test j2eq.L == 1.0
    end

    @testset "J2EqOE - Constructor from Vector" begin
        vec = [0.001, 0.1, 0.2, 0.05, 0.06, 1.0]
        j2eq = J2EqOE(vec)
        @test j2eq isa J2EqOE{Float64}
        @test j2eq.n == 0.001
        @test j2eq.L == 1.0
    end

    @testset "J2EqOE - Type Promotion" begin
        j2eq = J2EqOE(0.001, 0.1f0, 0.2, 0, 0.06, 1.0)
        @test j2eq isa J2EqOE
        @test eltype(j2eq) <: AbstractFloat
    end

    @testset "J2EqOE - params() Method" begin
        j2eq = J2EqOE(0.001, 0.1, 0.2, 0.05, 0.06, 1.0)
        p = params(j2eq)

        @test p isa SVector{6,Float64}
        @test p[1] == 0.001
        @test p[6] == 1.0
    end

    @testset "J2EqOE - Base.one()" begin
        j2eq_one = one(J2EqOE)
        @test j2eq_one isa J2EqOE{Float64}
        @test all(params(j2eq_one) .== 0.0)
    end

    @testset "J2EqOE - Base.getindex" begin
        j2eq = J2EqOE(0.001, 0.1, 0.2, 0.05, 0.06, 1.0)

        @test j2eq[1] == 0.001  # n
        @test j2eq[2] == 0.1    # h
        @test j2eq[3] == 0.2    # k
        @test j2eq[4] == 0.05   # p
        @test j2eq[5] == 0.06   # q
        @test j2eq[6] == 1.0    # L

        @test_throws BoundsError j2eq[0]
        @test_throws BoundsError j2eq[7]
    end

    @testset "Type Stability" begin
        μ = 398600.4418

        # ModEq Float64
        kep_f64 = Keplerian(7000.0, 0.1, 0.2, 0.3, 0.4, 0.5)
        modeq_f64 = ModEq(kep_f64, μ)
        @test eltype(modeq_f64) === Float64

        # ModEq Float32
        kep_f32 = Keplerian(7000.0f0, 0.1f0, 0.2f0, 0.3f0, 0.4f0, 0.5f0)
        modeq_f32 = ModEq(kep_f32, Float32(μ))
        @test eltype(modeq_f32) <: AbstractFloat
    end

    @testset "ModEqN Coordinate Set" begin
        @testset "Inner constructor ModEqN{T}" begin
            modeqn = ModEqN{Float64}(0.001, 0.1, 0.2, 0.05, 0.06, 1.0)
            @test modeqn isa ModEqN{Float64}
            @test modeqn.η == 0.001
            @test modeqn.f == 0.1
            @test modeqn.g == 0.2
            @test modeqn.h == 0.05
            @test modeqn.k == 0.06
            @test modeqn.L == 1.0
        end

        @testset "Constructor from vector" begin
            vec = [0.001, 0.1, 0.2, 0.05, 0.06, 1.0]
            modeqn = ModEqN(vec)
            @test modeqn isa ModEqN{Float64}
            @test modeqn.η == 0.001
            @test modeqn.L == 1.0

            # Test with Int vector
            vec_int = [1, 0, 0, 0, 0, 1]
            modeqn_int = ModEqN(vec_int)
            @test modeqn_int isa ModEqN{Int}
        end

        @testset "Type promotion" begin
            # Mixed types
            modeqn = ModEqN(0.001, 0.1f0, 0.2, 0, 0.06, 1.0)
            @test modeqn isa ModEqN
            @test eltype(modeqn) <: AbstractFloat

            # All Float32
            modeqn_f32 = ModEqN(0.001f0, 0.1f0, 0.2f0, 0.05f0, 0.06f0, 1.0f0)
            @test modeqn_f32 isa ModEqN{Float32}
        end

        @testset "StaticVector constructor" begin
            svec = SVector{6}(0.001, 0.1, 0.2, 0.05, 0.06, 1.0)
            modeqn = ModEqN(svec)
            @test modeqn isa ModEqN{Float64}
            @test modeqn.η == 0.001
        end

        @testset "Base.one()" begin
            modeqn_one = one(ModEqN)
            @test modeqn_one isa ModEqN{Float64}
            @test all(params(modeqn_one) .== 0.0)

            modeqn_one_f32 = one(ModEqN, T=Float32)
            @test modeqn_one_f32 isa ModEqN{Float32}
        end

        @testset "params() method" begin
            modeqn = ModEqN(0.001, 0.1, 0.2, 0.05, 0.06, 1.0)
            p = params(modeqn)

            @test p isa SVector{6,Float64}
            @test p[1] == 0.001
            @test p[2] == 0.1
            @test p[3] == 0.2
            @test p[4] == 0.05
            @test p[5] == 0.06
            @test p[6] == 1.0
        end

        @testset "Base.getindex" begin
            modeqn = ModEqN(0.001, 0.1, 0.2, 0.05, 0.06, 1.0)

            @test modeqn[1] == 0.001  # η
            @test modeqn[2] == 0.1    # f
            @test modeqn[3] == 0.2    # g
            @test modeqn[4] == 0.05   # h
            @test modeqn[5] == 0.06   # k
            @test modeqn[6] == 1.0    # L

            @test_throws BoundsError modeqn[0]
            @test_throws BoundsError modeqn[7]
        end
    end
end
