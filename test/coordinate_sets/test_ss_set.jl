using Test
using AstroCoords
using StaticArrays

@testset "StiefelScheifele Coordinate Set" begin
    @testset "Inner Constructor" begin
        # Test inner constructor (lines 36-37)
        α1, α2, α3, α4 = 1.0, 2.0, 3.0, 4.0
        β1, β2, β3, β4 = 0.1, 0.2, 0.3, 0.4
        ω = 0.5
        t = 100.0
        
        ss = StiefelScheifele{Float64}(α1, α2, α3, α4, β1, β2, β3, β4, ω, t)
        
        @test ss isa StiefelScheifele{Float64}
        @test ss.α1 == 1.0
        @test ss.α2 == 2.0
        @test ss.α3 == 3.0
        @test ss.α4 == 4.0
        @test ss.β1 == 0.1
        @test ss.β2 == 0.2
        @test ss.β3 == 0.3
        @test ss.β4 == 0.4
        @test ss.ω == 0.5
        @test ss.t == 100.0
    end
    
    @testset "Constructor from AbstractVector" begin
        # Test constructor from AbstractVector (lines 45-46)
        vec = [1.0, 2.0, 3.0, 4.0, 0.1, 0.2, 0.3, 0.4, 0.5, 100.0]
        ss = StiefelScheifele(vec)
        
        @test ss isa StiefelScheifele{Float64}
        @test ss.α1 == 1.0
        @test ss.α2 == 2.0
        @test ss.α3 == 3.0
        @test ss.α4 == 4.0
        @test ss.β1 == 0.1
        @test ss.β2 == 0.2
        @test ss.β3 == 0.3
        @test ss.β4 == 0.4
        @test ss.ω == 0.5
        @test ss.t == 100.0
        
        # Test with SVector
        svec = SVector{10}(1.0, 2.0, 3.0, 4.0, 0.1, 0.2, 0.3, 0.4, 0.5, 100.0)
        ss2 = StiefelScheifele(svec)
        @test ss2 isa StiefelScheifele{Float64}
        @test ss2.α1 == 1.0
    end
    
    @testset "Constructor with Individual Parameters" begin
        # Test constructor with individual parameters (line 48)
        ss = StiefelScheifele(1.0, 2.0, 3.0, 4.0, 0.1, 0.2, 0.3, 0.4, 0.5, 100.0)
        
        @test ss isa StiefelScheifele{Float64}
        @test ss.α1 == 1.0
        @test ss.α2 == 2.0
        @test ss.α3 == 3.0
        @test ss.α4 == 4.0
        @test ss.β1 == 0.1
        @test ss.β2 == 0.2
        @test ss.β3 == 0.3
        @test ss.β4 == 0.4
        @test ss.ω == 0.5
        @test ss.t == 100.0
    end
    
    @testset "Constructor with Type Promotion" begin
        # Test type promotion (lines 51-52)
        α1, α2, α3, α4 = 1, 2.0, 3, 4.0  # Mix Int and Float64
        β1, β2, β3, β4 = 0.1f0, 0.2, 0.3f0, 0.4  # Mix Float32 and Float64
        ω = 0.5
        t = 100
        
        ss = StiefelScheifele(α1, α2, α3, α4, β1, β2, β3, β4, ω, t)
        
        # Should promote to Float64 (highest precision)
        @test ss isa StiefelScheifele{Float64}
        @test typeof(ss.α1) == Float64
        @test typeof(ss.β1) == Float64
        @test typeof(ss.ω) == Float64
        @test typeof(ss.t) == Float64
    end
    
    @testset "Element Access - All 10 Elements" begin
        # Test all element access (lines 59-60)
        ss = StiefelScheifele(1.0, 2.0, 3.0, 4.0, 0.1, 0.2, 0.3, 0.4, 0.5, 100.0)
        
        # Test indexed access
        @test ss[1] == 1.0
        @test ss[2] == 2.0
        @test ss[3] == 3.0
        @test ss[4] == 4.0
        @test ss[5] == 0.1
        @test ss[6] == 0.2
        @test ss[7] == 0.3
        @test ss[8] == 0.4
        @test ss[9] == 0.5
        @test ss[10] == 100.0
        
        # Test out of bounds
        @test_throws BoundsError ss[0]
        @test_throws BoundsError ss[11]
        @test_throws BoundsError ss[-1]
    end
    
    @testset "Struct Consistency with SS Transformations" begin
        # Verify struct fields match transformation expectations
        μ = 398600.4418  # Earth μ (km³/s²)
        cart = Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)
        
        # Create config for SS transformation
        config = RegularizedCoordinateConfig(params(cart), μ)
        ϕ = compute_initial_phi(params(cart), μ, config)
        
        # Transform to SS coordinates
        ss_vec = cart2StiefelScheifele(params(cart), μ, ϕ, config)
        ss = StiefelScheifele(ss_vec)
        
        # Verify all 10 elements are accessible and match
        @test ss.α1 == ss_vec[1]
        @test ss.α2 == ss_vec[2]
        @test ss.α3 == ss_vec[3]
        @test ss.α4 == ss_vec[4]
        @test ss.β1 == ss_vec[5]
        @test ss.β2 == ss_vec[6]
        @test ss.β3 == ss_vec[7]
        @test ss.β4 == ss_vec[8]
        @test ss.ω == ss_vec[9]
        @test ss.t == ss_vec[10]
        
        # Verify params() returns all elements
        p = params(ss)
        @test length(p) == 10
        @test p[1] == ss.α1
        @test p[10] == ss.t
    end
    
    @testset "Various α and β Values" begin
        # Test with different value patterns
        
        # All zeros
        ss_zero = StiefelScheifele(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        @test ss_zero.α1 == 0.0
        @test ss_zero.ω == 0.0
        
        # All ones
        ss_one = StiefelScheifele(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
        @test ss_one.α1 == 1.0
        @test ss_one.β1 == 1.0
        
        # Negative values
        ss_neg = StiefelScheifele(-1.0, -2.0, -3.0, -4.0, -0.1, -0.2, -0.3, -0.4, -0.5, -100.0)
        @test ss_neg.α1 == -1.0
        @test ss_neg.ω == -0.5
        @test ss_neg.t == -100.0
        
        # Large values
        ss_large = StiefelScheifele(1e6, 1e6, 1e6, 1e6, 1e3, 1e3, 1e3, 1e3, 1e6, 1e9)
        @test ss_large.α1 == 1e6
        @test ss_large.ω == 1e6
        
        # Small values
        ss_small = StiefelScheifele(1e-10, 1e-10, 1e-10, 1e-10, 1e-12, 1e-12, 1e-12, 1e-12, 1e-15, 1e-6)
        @test ss_small.α1 == 1e-10
        @test ss_small.t == 1e-6
    end
    
    @testset "params() Function" begin
        # Test params extraction
        ss = StiefelScheifele(1.0, 2.0, 3.0, 4.0, 0.1, 0.2, 0.3, 0.4, 0.5, 100.0)
        p = params(ss)
        
        @test p isa SVector{10, Float64}
        @test p[1] == ss.α1
        @test p[2] == ss.α2
        @test p[3] == ss.α3
        @test p[4] == ss.α4
        @test p[5] == ss.β1
        @test p[6] == ss.β2
        @test p[7] == ss.β3
        @test p[8] == ss.β4
        @test p[9] == ss.ω
        @test p[10] == ss.t
    end
    
    @testset "Base.one Constructor" begin
        # Test zero initialization
        ss_one = one(StiefelScheifele)
        @test ss_one isa StiefelScheifele{Float64}
        @test ss_one.α1 == 0.0
        @test ss_one.α2 == 0.0
        @test ss_one.ω == 0.0
        @test ss_one.t == 0.0
        
        # Test with specific type
        ss_one_f32 = one(StiefelScheifele, T=Float32)
        @test ss_one_f32 isa StiefelScheifele{Float32}
        @test typeof(ss_one_f32.α1) == Float32
    end
end
