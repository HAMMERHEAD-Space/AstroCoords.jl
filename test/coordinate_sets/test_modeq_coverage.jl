using Test
using AstroCoords
using StaticArrays

@testset "ModEq Constructor Coverage" begin
    # Target: Lines 30, 35-37, 40-41, 44, 47-48
    # These are constructor variants - no literature needed (API tests)
    
    @testset "AbstractVector constructor" begin
        # Line 30: ModEq(X::AbstractVector{T})
        vec = [7000.0, 0.01, 0.02, 0.001, 0.002, 0.5]
        meq = ModEq(vec)
        @test meq isa ModEq{Float64}
        @test meq.p ≈ 7000.0
        @test meq.f ≈ 0.01
        @test meq.g ≈ 0.02
        @test meq.h ≈ 0.001
        @test meq.k ≈ 0.002
        @test meq.L ≈ 0.5
    end
    
    @testset "Individual args with type promotion" begin
        # Lines 35-37: ModEq(p::P, f::F, g::G, h::H, k::K, L::LT) with promote_type
        # Test Int → Float64 promotion
        meq1 = ModEq(7000, 0.01, 0.02, 0.001, 0.002, 0.5)
        @test meq1 isa ModEq{Float64}
        @test meq1.p ≈ 7000.0
        
        # Test Float32 → Float64 promotion
        meq2 = ModEq(7000.0f0, 0.01f0, 0.02, 0.001, 0.002, 0.5)
        @test meq2 isa ModEq{Float64}
        @test meq2.p ≈ 7000.0
        
        # Test mixed Int/Float32/Float64
        meq3 = ModEq(7000, 0.01f0, 0.02, 0.001f0, 0.002, 0.5f0)
        @test meq3 isa ModEq{Float64}
        @test meq3.f ≈ 0.01
    end
    
    @testset "StaticVector constructor" begin
        # Lines 40-41: ModEq(g::StaticVector{N,T})
        sv = SVector{6,Float64}(7000.0, 0.01, 0.02, 0.001, 0.002, 0.5)
        meq = ModEq(sv)
        @test meq isa ModEq{Float64}
        @test meq.p ≈ 7000.0
        @test meq.L ≈ 0.5
        
        # Test with Float32 StaticVector
        sv32 = SVector{6,Float32}(7000.0f0, 0.01f0, 0.02f0, 0.001f0, 0.002f0, 0.5f0)
        meq32 = ModEq(sv32)
        @test meq32 isa ModEq{Float32}
        @test meq32.p ≈ 7000.0f0
    end
    
    @testset "ModEq{T}(StaticVector) type conversion" begin
        # Line 41: ModEq{T}(g::StaticVector)
        sv = SVector{6,Float32}(7000.0f0, 0.01f0, 0.02f0, 0.001f0, 0.002f0, 0.5f0)
        meq = ModEq{Float64}(sv)
        @test meq isa ModEq{Float64}
        @test meq.p ≈ 7000.0
    end
    
    @testset "Base.one() with default type" begin
        # Lines 47-48: Base.one(::Type{M}; T::DataType=Float64)
        meq_one = Base.one(ModEq)
        @test meq_one isa ModEq{Float64}
        @test meq_one.p == 0.0
        @test meq_one.f == 0.0
        @test meq_one.g == 0.0
        @test meq_one.h == 0.0
        @test meq_one.k == 0.0
        @test meq_one.L == 0.0
    end
    
    @testset "Base.one() with custom type" begin
        # Lines 47-48: Base.one with T=Float32
        meq_one32 = Base.one(ModEq; T=Float32)
        @test meq_one32 isa ModEq{Float32}
        @test meq_one32.p == 0.0f0
        @test meq_one32.f == 0.0f0
    end
    
    @testset "Field access via indexing" begin
        # Lines 51-69: Base.getindex(p::ModEq{T}, i::Int)
        meq = ModEq(7000.0, 0.01, 0.02, 0.001, 0.002, 0.5)
        @test meq[1] ≈ 7000.0  # p
        @test meq[2] ≈ 0.01    # f
        @test meq[3] ≈ 0.02    # g
        @test meq[4] ≈ 0.001   # h
        @test meq[5] ≈ 0.002   # k
        @test meq[6] ≈ 0.5     # L
        @test_throws BoundsError meq[7]
        @test_throws BoundsError meq[0]
    end
    
    @testset "params() extraction" begin
        # Line 44: params(g::ModEq{T})
        meq = ModEq(7000.0, 0.01, 0.02, 0.001, 0.002, 0.5)
        p = params(meq)
        @test p isa SVector{6,Float64}
        @test p[1] ≈ 7000.0
        @test p[6] ≈ 0.5
    end
    
    @testset "Near-circular orbit edge case" begin
        # Test near-circular orbit: f,g,h,k ≈ 0
        meq_circ = ModEq(7000.0, 1e-10, 1e-10, 1e-12, 1e-12, 0.5)
        @test meq_circ.f ≈ 1e-10 atol=1e-15
        @test meq_circ.g ≈ 1e-10 atol=1e-15
        @test meq_circ.h ≈ 1e-12 atol=1e-15
        @test meq_circ.k ≈ 1e-12 atol=1e-15
    end
    
    @testset "Equatorial orbit edge case" begin
        # Test equatorial orbit: h,k ≈ 0
        meq_equat = ModEq(7000.0, 0.01, 0.02, 1e-12, 1e-12, 0.5)
        @test meq_equat.h ≈ 1e-12 atol=1e-15
        @test meq_equat.k ≈ 1e-12 atol=1e-15
        @test meq_equat.f ≈ 0.01
        @test meq_equat.g ≈ 0.02
    end
end
