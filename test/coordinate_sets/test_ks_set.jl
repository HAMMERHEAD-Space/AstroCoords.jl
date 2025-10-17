using Test
using AstroCoords
using StaticArrays

@testset "KustaanheimoStiefel Coordinate Set" begin
    @testset "Inner Constructor" begin
        # Test inner constructor (lines 37-40)
        u₁, u₂, u₃, u₄ = 1.0, 2.0, 3.0, 4.0
        u₁_prime, u₂_prime, u₃_prime, u₄_prime = 0.1, 0.2, 0.3, 0.4
        h = -0.5
        τ = 100.0
        
        ks = KustaanheimoStiefel{Float64}(u₁, u₂, u₃, u₄, u₁_prime, u₂_prime, u₃_prime, u₄_prime, h, τ)
        
        @test ks isa KustaanheimoStiefel{Float64}
        @test ks.u₁ == 1.0
        @test ks.u₂ == 2.0
        @test ks.u₃ == 3.0
        @test ks.u₄ == 4.0
        @test ks.u₁_prime == 0.1
        @test ks.u₂_prime == 0.2
        @test ks.u₃_prime == 0.3
        @test ks.u₄_prime == 0.4
        @test ks.h == -0.5
        @test ks.τ == 100.0
    end
    
    @testset "Constructor from AbstractVector" begin
        # Test constructor from AbstractVector (lines 48-49)
        vec = [1.0, 2.0, 3.0, 4.0, 0.1, 0.2, 0.3, 0.4, -0.5, 100.0]
        ks = KustaanheimoStiefel(vec)
        
        @test ks isa KustaanheimoStiefel{Float64}
        @test ks.u₁ == 1.0
        @test ks.u₂ == 2.0
        @test ks.u₃ == 3.0
        @test ks.u₄ == 4.0
        @test ks.u₁_prime == 0.1
        @test ks.u₂_prime == 0.2
        @test ks.u₃_prime == 0.3
        @test ks.u₄_prime == 0.4
        @test ks.h == -0.5
        @test ks.τ == 100.0
        
        # Test with SVector
        svec = SVector{10}(1.0, 2.0, 3.0, 4.0, 0.1, 0.2, 0.3, 0.4, -0.5, 100.0)
        ks2 = KustaanheimoStiefel(svec)
        @test ks2 isa KustaanheimoStiefel{Float64}
        @test ks2.u₁ == 1.0
    end
    
    @testset "Constructor with Type Promotion" begin
        # Test type promotion (line 52)
        u₁, u₂, u₃, u₄ = 1, 2.0, 3, 4.0  # Mix Int and Float64
        u₁_prime, u₂_prime, u₃_prime, u₄_prime = 0.1f0, 0.2, 0.3f0, 0.4  # Mix Float32 and Float64
        h = -0.5
        τ = 100
        
        ks = KustaanheimoStiefel(u₁, u₂, u₃, u₄, u₁_prime, u₂_prime, u₃_prime, u₄_prime, h, τ)
        
        # Should promote to Float64 (highest precision)
        @test ks isa KustaanheimoStiefel{Float64}
        @test typeof(ks.u₁) == Float64
        @test typeof(ks.u₁_prime) == Float64
        @test typeof(ks.h) == Float64
        @test typeof(ks.τ) == Float64
    end
    
    @testset "Element Access - All 10 Elements" begin
        # Test all element access (lines 64-75)
        ks = KustaanheimoStiefel(1.0, 2.0, 3.0, 4.0, 0.1, 0.2, 0.3, 0.4, -0.5, 100.0)
        
        # Test indexed access (lines 64-75)
        @test ks[1] == 1.0
        @test ks[2] == 2.0
        @test ks[3] == 3.0
        @test ks[4] == 4.0
        @test ks[5] == 0.1
        @test ks[6] == 0.2
        @test ks[7] == 0.3
        @test ks[8] == 0.4
        @test ks[9] == -0.5
        @test ks[10] == 100.0
        
        # Test out of bounds
        @test_throws BoundsError ks[0]
        @test_throws BoundsError ks[11]
        @test_throws BoundsError ks[-1]
    end
    
    @testset "Struct Consistency with KS Transformations" begin
        # Verify struct fields match transformation expectations
        μ = 398600.4418  # Earth μ (km³/s²)
        cart = Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)
        
        # Create config for KS transformation
        config = RegularizedCoordinateConfig(params(cart), μ)
        
        # Transform to KS coordinates
        ks_vec = cart2KS(params(cart), μ, config)
        ks = KustaanheimoStiefel(ks_vec)
        
        # Verify all 10 elements are accessible and match
        @test ks.u₁ == ks_vec[1]
        @test ks.u₂ == ks_vec[2]
        @test ks.u₃ == ks_vec[3]
        @test ks.u₄ == ks_vec[4]
        @test ks.u₁_prime == ks_vec[5]
        @test ks.u₂_prime == ks_vec[6]
        @test ks.u₃_prime == ks_vec[7]
        @test ks.u₄_prime == ks_vec[8]
        @test ks.h == ks_vec[9]
        @test ks.τ == ks_vec[10]
        
        # Verify params() returns all elements
        p = params(ks)
        @test length(p) == 10
        @test p[1] == ks.u₁
        @test p[10] == ks.τ
    end
    
    @testset "Various u and u_prime Values" begin
        # Test with different value patterns
        
        # All zeros
        ks_zero = KustaanheimoStiefel(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        @test ks_zero.u₁ == 0.0
        @test ks_zero.h == 0.0
        
        # All ones
        ks_one = KustaanheimoStiefel(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        @test ks_one.u₁ == 1.0
        @test ks_one.u₁_prime == 1.0
        
        # Negative values
        ks_neg = KustaanheimoStiefel(-1.0, -2.0, -3.0, -4.0, -0.1, -0.2, -0.3, -0.4, -0.5, -100.0)
        @test ks_neg.u₁ == -1.0
        @test ks_neg.h == -0.5
        @test ks_neg.τ == -100.0
        
        # Large values
        ks_large = KustaanheimoStiefel(1e6, 1e6, 1e6, 1e6, 1e3, 1e3, 1e3, 1e3, -1e6, 1e9)
        @test ks_large.u₁ == 1e6
        @test ks_large.h == -1e6
        
        # Small values
        ks_small = KustaanheimoStiefel(1e-10, 1e-10, 1e-10, 1e-10, 1e-12, 1e-12, 1e-12, 1e-12, -1e-15, 1e-6)
        @test ks_small.u₁ == 1e-10
        @test ks_small.τ == 1e-6
    end
    
    @testset "params() Function" begin
        # Test params extraction
        ks = KustaanheimoStiefel(1.0, 2.0, 3.0, 4.0, 0.1, 0.2, 0.3, 0.4, -0.5, 100.0)
        p = params(ks)
        
        @test p isa SVector{10, Float64}
        @test p[1] == ks.u₁
        @test p[2] == ks.u₂
        @test p[3] == ks.u₃
        @test p[4] == ks.u₄
        @test p[5] == ks.u₁_prime
        @test p[6] == ks.u₂_prime
        @test p[7] == ks.u₃_prime
        @test p[8] == ks.u₄_prime
        @test p[9] == ks.h
        @test p[10] == ks.τ
    end
    
    @testset "Base.one Constructor" begin
        # Test zero initialization
        ks_one = one(KustaanheimoStiefel)
        @test ks_one isa KustaanheimoStiefel{Float64}
        @test ks_one.u₁ == 0.0
        @test ks_one.u₂ == 0.0
        @test ks_one.h == 0.0
        @test ks_one.τ == 0.0
        
        # Test with specific type
        ks_one_f32 = one(KustaanheimoStiefel, T=Float32)
        @test ks_one_f32 isa KustaanheimoStiefel{Float32}
        @test typeof(ks_one_f32.u₁) == Float32
    end
end
