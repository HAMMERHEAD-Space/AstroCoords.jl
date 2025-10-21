using Test
using AstroCoords
using StaticArrays

@testset "Keplerian Constructor Coverage" begin
    # Target: Lines 31, 36-38, 41-42, 44, 47, 50-51
    # These are constructor variants - no literature needed (API tests)
    
    @testset "AbstractVector constructor" begin
        # Line 31: Keplerian(X::AbstractVector{T})
        vec = [7000.0, 0.01, 0.1, 0.2, 0.3, 0.5]
        kep = Keplerian(vec)
        @test kep isa Keplerian{Float64}
        @test kep.a ≈ 7000.0
        @test kep.e ≈ 0.01
        @test kep.i ≈ 0.1
        @test kep.Ω ≈ 0.2
        @test kep.ω ≈ 0.3
        @test kep.f ≈ 0.5
    end
    
    @testset "Individual args with type promotion" begin
        # Lines 36-38: Keplerian(a::A, e::E, i::I, Ω::O, ω::W, f::F) with promote_type
        # Test Int → Float64 promotion
        kep1 = Keplerian(7000, 0.01, 0.1, 0.2, 0.3, 0.5)
        @test kep1 isa Keplerian{Float64}
        @test kep1.a ≈ 7000.0
        
        # Test Float32 → Float64 promotion
        kep2 = Keplerian(7000.0f0, 0.01f0, 0.1, 0.2, 0.3, 0.5)
        @test kep2 isa Keplerian{Float64}
        @test kep2.a ≈ 7000.0
        
        # Test mixed Int/Float32/Float64
        kep3 = Keplerian(7000, 0.01f0, 0.1, 0.2f0, 0.3, 0.5f0)
        @test kep3 isa Keplerian{Float64}
        @test kep3.e ≈ 0.01
    end
    
    @testset "StaticVector{6} constructor" begin
        # Lines 41-42, 44: Keplerian(g::StaticVector{N,T})
        sv = SVector{6,Float64}(7000.0, 0.01, 0.1, 0.2, 0.3, 0.5)
        kep = Keplerian(sv)
        @test kep isa Keplerian{Float64}
        @test kep.a ≈ 7000.0
        @test kep.f ≈ 0.5
        
        # Test with Float32 StaticVector
        sv32 = SVector{6,Float32}(7000.0f0, 0.01f0, 0.1f0, 0.2f0, 0.3f0, 0.5f0)
        kep32 = Keplerian(sv32)
        @test kep32 isa Keplerian{Float32}
        @test kep32.a ≈ 7000.0f0
    end
    
    @testset "Keplerian{T}(StaticVector) type conversion" begin
        # Line 44: Keplerian{T}(g::StaticVector)
        sv = SVector{6,Float32}(7000.0f0, 0.01f0, 0.1f0, 0.2f0, 0.3f0, 0.5f0)
        kep = Keplerian{Float64}(sv)
        @test kep isa Keplerian{Float64}
        @test kep.a ≈ 7000.0
    end
    
    @testset "Base.one() with default type" begin
        # Lines 47, 50-51: Base.one(::Type{K}; T::DataType=Float64)
        kep_one = Base.one(Keplerian)
        @test kep_one isa Keplerian{Float64}
        @test kep_one.a == 0.0
        @test kep_one.e == 0.0
        @test kep_one.i == 0.0
        @test kep_one.Ω == 0.0
        @test kep_one.ω == 0.0
        @test kep_one.f == 0.0
    end
    
    @testset "Base.one() with custom type" begin
        # Lines 47, 50-51: Base.one with T=Float32
        kep_one32 = Base.one(Keplerian; T=Float32)
        @test kep_one32 isa Keplerian{Float32}
        @test kep_one32.a == 0.0f0
        @test kep_one32.e == 0.0f0
    end
    
    @testset "Field access via indexing" begin
        # Lines 54-68: Base.getindex(p::Keplerian{T}, i::Int)
        kep = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.5)
        @test kep[1] ≈ 7000.0  # a
        @test kep[2] ≈ 0.01    # e
        @test kep[3] ≈ 0.1     # i
        @test kep[4] ≈ 0.2     # Ω
        @test kep[5] ≈ 0.3     # ω
        @test kep[6] ≈ 0.5     # f
        @test_throws BoundsError kep[7]
        @test_throws BoundsError kep[0]
    end
    
    @testset "params() extraction" begin
        # Line 46: params(g::Keplerian{T})
        kep = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.5)
        p = params(kep)
        @test p isa SVector{6,Float64}
        @test p[1] ≈ 7000.0
        @test p[6] ≈ 0.5
    end
end
