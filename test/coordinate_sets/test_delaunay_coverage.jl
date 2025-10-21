using Test
using AstroCoords
using StaticArrays

@testset "Delaunay Constructor Coverage" begin
    # Target: Lines 31, 36-38, 41-42, 45, 48-49
    # These are constructor variants - no literature needed (API tests)
    
    @testset "AbstractVector constructor" begin
        # Line 31: Delaunay(X::AbstractVector{T})
        vec = [1000.0, 950.0, 900.0, 0.5, 0.3, 0.2]
        del = Delaunay(vec)
        @test del isa Delaunay{Float64}
        @test del.L ≈ 1000.0
        @test del.G ≈ 950.0
        @test del.H ≈ 900.0
        @test del.M ≈ 0.5
        @test del.ω ≈ 0.3
        @test del.Ω ≈ 0.2
    end
    
    @testset "Individual args with type promotion" begin
        # Lines 36-38: Delaunay(L::LT, G::GT, H::HT, M::MT, ω::PT, Ω::OmT) with promote_type
        # Test Int → Float64 promotion
        del1 = Delaunay(1000, 950, 900, 0.5, 0.3, 0.2)
        @test del1 isa Delaunay{Float64}
        @test del1.L ≈ 1000.0
        
        # Test Float32 → Float64 promotion
        del2 = Delaunay(1000.0f0, 950.0f0, 900.0, 0.5, 0.3, 0.2)
        @test del2 isa Delaunay{Float64}
        @test del2.L ≈ 1000.0
        
        # Test mixed Int/Float32/Float64
        del3 = Delaunay(1000, 950.0f0, 900, 0.5f0, 0.3, 0.2f0)
        @test del3 isa Delaunay{Float64}
        @test del3.G ≈ 950.0
    end
    
    @testset "StaticVector constructor" begin
        # Lines 41-42: Delaunay(g::StaticVector{N,T})
        sv = SVector{6,Float64}(1000.0, 950.0, 900.0, 0.5, 0.3, 0.2)
        del = Delaunay(sv)
        @test del isa Delaunay{Float64}
        @test del.L ≈ 1000.0
        @test del.Ω ≈ 0.2
        
        # Test with Float32 StaticVector
        sv32 = SVector{6,Float32}(1000.0f0, 950.0f0, 900.0f0, 0.5f0, 0.3f0, 0.2f0)
        del32 = Delaunay(sv32)
        @test del32 isa Delaunay{Float32}
        @test del32.L ≈ 1000.0f0
    end
    
    @testset "Delaunay{T}(StaticVector) type conversion" begin
        # Line 42: Delaunay{T}(g::StaticVector)
        sv = SVector{6,Float32}(1000.0f0, 950.0f0, 900.0f0, 0.5f0, 0.3f0, 0.2f0)
        del = Delaunay{Float64}(sv)
        @test del isa Delaunay{Float64}
        @test del.L ≈ 1000.0
    end
    
    @testset "Base.one() with default type" begin
        # Lines 48-49: Base.one(::Type{D}; T::DataType=Float64)
        del_one = Base.one(Delaunay)
        @test del_one isa Delaunay{Float64}
        @test del_one.L == 0.0
        @test del_one.G == 0.0
        @test del_one.H == 0.0
        @test del_one.M == 0.0
        @test del_one.ω == 0.0
        @test del_one.Ω == 0.0
    end
    
    @testset "Base.one() with custom type" begin
        # Lines 48-49: Base.one with T=Float32
        del_one32 = Base.one(Delaunay; T=Float32)
        @test del_one32 isa Delaunay{Float32}
        @test del_one32.L == 0.0f0
        @test del_one32.G == 0.0f0
    end
    
    @testset "Field access via indexing" begin
        # Lines 52-69: Base.getindex(p::Delaunay{T}, i::Int)
        del = Delaunay(1000.0, 950.0, 900.0, 0.5, 0.3, 0.2)
        @test del[1] ≈ 1000.0  # L
        @test del[2] ≈ 950.0   # G
        @test del[3] ≈ 900.0   # H
        @test del[4] ≈ 0.5     # M
        @test del[5] ≈ 0.3     # ω
        @test del[6] ≈ 0.2     # Ω
        @test_throws BoundsError del[7]
        @test_throws BoundsError del[0]
    end
    
    @testset "params() extraction" begin
        # Line 45: params(g::Delaunay{T})
        del = Delaunay(1000.0, 950.0, 900.0, 0.5, 0.3, 0.2)
        p = params(del)
        @test p isa SVector{6,Float64}
        @test p[1] ≈ 1000.0
        @test p[6] ≈ 0.2
    end
    
    @testset "Circular orbit edge case" begin
        # Test circular orbit: G = L (eccentricity = 0)
        del_circ = Delaunay(1000.0, 1000.0, 900.0, 0.5, 0.3, 0.2)
        @test del_circ.L ≈ del_circ.G
        # G/L = 1.0 → e = 0
        # Note: This tests constructor handles G=L without issues
    end
    
    @testset "Equatorial orbit edge case" begin
        # Test equatorial orbit: H = G (inclination = 0)
        del_equat = Delaunay(1000.0, 950.0, 950.0, 0.5, 0.3, 0.2)
        @test del_equat.G ≈ del_equat.H
        # H = G → cos(i) = 1 → i = 0
    end
    
    @testset "Type promotion with mixed numeric types" begin
        # Test complex type promotion scenario
        del_mixed = Delaunay(Int64(1000), Float32(950.0), Float64(900.0), 
                              Float32(0.5), Int32(0), Float64(0.2))
        @test del_mixed isa Delaunay{Float64}
        @test del_mixed.L ≈ 1000.0
        @test del_mixed.G ≈ 950.0
        @test del_mixed.H ≈ 900.0
        @test del_mixed.M ≈ 0.5
        @test del_mixed.ω ≈ 0.0
        @test del_mixed.Ω ≈ 0.2
    end
end
