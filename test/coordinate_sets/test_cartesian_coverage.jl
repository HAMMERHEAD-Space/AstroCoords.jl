using Test
using AstroCoords
using StaticArrays

@testset "Cartesian Constructor Coverage" begin
    # Target: Lines 26, 33-35, 38-39, 41, 44, 47-48
    # These are constructor variants - no literature needed (API tests)
    
    @testset "AbstractVector constructor" begin
        # Line 26: Cartesian(X::AbstractVector{T})
        vec = [7000.0, 1000.0, 500.0, 1.5, 2.0, 3.5]
        cart = Cartesian(vec)
        @test cart isa Cartesian{Float64}
        @test cart.x ≈ 7000.0
        @test cart.y ≈ 1000.0
        @test cart.z ≈ 500.0
        @test cart.ẋ ≈ 1.5
        @test cart.ẏ ≈ 2.0
        @test cart.ż ≈ 3.5
    end
    
    @testset "Individual args with type promotion" begin
        # Lines 33-35: Cartesian(x::X, y::Y, z::Z, ẋ::XV, ẏ::YV, ż::ZV) with promote_type
        # Test Int → Float64 promotion
        cart1 = Cartesian(7000, 1000, 500, 1.5, 2.0, 3.5)
        @test cart1 isa Cartesian{Float64}
        @test cart1.x ≈ 7000.0
        
        # Test Float32 → Float64 promotion
        cart2 = Cartesian(7000.0f0, 1000.0f0, 500.0, 1.5, 2.0, 3.5)
        @test cart2 isa Cartesian{Float64}
        @test cart2.x ≈ 7000.0
        
        # Test mixed Int/Float32/Float64
        cart3 = Cartesian(7000, 1000.0f0, 500, 1.5f0, 2.0, 3.5f0)
        @test cart3 isa Cartesian{Float64}
        @test cart3.y ≈ 1000.0
    end
    
    @testset "StaticVector constructor" begin
        # Lines 38-39, 41: Cartesian(g::StaticVector{N,T})
        sv = SVector{6,Float64}(7000.0, 1000.0, 500.0, 1.5, 2.0, 3.5)
        cart = Cartesian(sv)
        @test cart isa Cartesian{Float64}
        @test cart.x ≈ 7000.0
        @test cart.ż ≈ 3.5
        
        # Test with Float32 StaticVector
        sv32 = SVector{6,Float32}(7000.0f0, 1000.0f0, 500.0f0, 1.5f0, 2.0f0, 3.5f0)
        cart32 = Cartesian(sv32)
        @test cart32 isa Cartesian{Float32}
        @test cart32.x ≈ 7000.0f0
    end
    
    @testset "Cartesian{T}(StaticVector) type conversion" begin
        # Line 41: Cartesian{T}(g::StaticVector)
        sv = SVector{6,Float32}(7000.0f0, 1000.0f0, 500.0f0, 1.5f0, 2.0f0, 3.5f0)
        cart = Cartesian{Float64}(sv)
        @test cart isa Cartesian{Float64}
        @test cart.x ≈ 7000.0
    end
    
    @testset "Base.one() with default type" begin
        # Lines 44, 47-48: Base.one(::Type{C}; T::DataType=Float64)
        cart_one = Base.one(Cartesian)
        @test cart_one isa Cartesian{Float64}
        @test cart_one.x == 0.0
        @test cart_one.y == 0.0
        @test cart_one.z == 0.0
        @test cart_one.ẋ == 0.0
        @test cart_one.ẏ == 0.0
        @test cart_one.ż == 0.0
    end
    
    @testset "Base.one() with custom type" begin
        # Lines 44, 47-48: Base.one with T=Float32
        cart_one32 = Base.one(Cartesian; T=Float32)
        @test cart_one32 isa Cartesian{Float32}
        @test cart_one32.x == 0.0f0
        @test cart_one32.ẋ == 0.0f0
    end
    
    @testset "Field access via indexing" begin
        # Lines 51-69: Base.getindex(p::Cartesian{T}, i::Int)
        cart = Cartesian(7000.0, 1000.0, 500.0, 1.5, 2.0, 3.5)
        @test cart[1] ≈ 7000.0  # x
        @test cart[2] ≈ 1000.0  # y
        @test cart[3] ≈ 500.0   # z
        @test cart[4] ≈ 1.5     # ẋ
        @test cart[5] ≈ 2.0     # ẏ
        @test cart[6] ≈ 3.5     # ż
        @test_throws BoundsError cart[7]
        @test_throws BoundsError cart[0]
    end
    
    @testset "params() extraction" begin
        # Line 43: params(g::Cartesian{T})
        cart = Cartesian(7000.0, 1000.0, 500.0, 1.5, 2.0, 3.5)
        p = params(cart)
        @test p isa SVector{6,Float64}
        @test p[1] ≈ 7000.0
        @test p[6] ≈ 3.5
    end
end
