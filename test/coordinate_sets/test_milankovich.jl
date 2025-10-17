using Test
using AstroCoords
using StaticArrays

@testset "Milankovich Coordinate Set" begin
    @testset "Inner Constructor" begin
        # Test inner constructor with explicit type
        mil = Milankovich{Float64}(1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 1.5)
        @test mil isa Milankovich{Float64}
        @test mil.hx == 1.0
        @test mil.hy == 2.0
        @test mil.hz == 3.0
        @test mil.ex == 0.1
        @test mil.ey == 0.2
        @test mil.ez == 0.3
        @test mil.L == 1.5
    end
    
    @testset "Constructor from AbstractVector" begin
        # Test constructor from vector
        vec = [1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 1.5]
        mil = Milankovich(vec)
        @test mil isa Milankovich{Float64}
        @test mil.hx == 1.0
        @test mil.hy == 2.0
        @test mil.hz == 3.0
        @test mil.ex == 0.1
        @test mil.ey == 0.2
        @test mil.ez == 0.3
        @test mil.L == 1.5
        
        # Test with different element types
        vec_int = [1, 2, 3, 0, 0, 0, 1]
        mil_int = Milankovich(vec_int)
        @test mil_int isa Milankovich{Int}
    end
    
    @testset "Constructor with Type Promotion" begin
        # Test type promotion with mixed types
        mil = Milankovich(1.0, 2, 3.0f0, 0.1, 0.2, 0.3, 1.5)
        @test mil isa Milankovich
        @test eltype(mil) <: AbstractFloat
        
        # Test with all Float32
        mil_f32 = Milankovich(1.0f0, 2.0f0, 3.0f0, 0.1f0, 0.2f0, 0.3f0, 1.5f0)
        @test mil_f32 isa Milankovich{Float32}
        
        # Test with mixed Int and Float
        mil_mixed = Milankovich(1, 2.0, 3, 0.1, 0, 0.3, 1)
        @test eltype(mil_mixed) <: AbstractFloat
    end
    
    @testset "StaticVector Constructor" begin
        # Test constructor from StaticVector
        svec = SVector{7}(1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 1.5)
        mil = Milankovich(svec)
        @test mil isa Milankovich{Float64}
        @test mil.hx == 1.0
        @test mil.hy == 2.0
        @test mil.hz == 3.0
        @test mil.ex == 0.1
        @test mil.ey == 0.2
        @test mil.ez == 0.3
        @test mil.L == 1.5
    end
    
    @testset "params() Method" begin
        # Test params() returns correct SVector
        mil = Milankovich(1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 1.5)
        p = params(mil)
        
        @test p isa SVector{7, Float64}
        @test p[1] == 1.0
        @test p[2] == 2.0
        @test p[3] == 3.0
        @test p[4] == 0.1
        @test p[5] == 0.2
        @test p[6] == 0.3
        @test p[7] == 1.5
    end
    
    @testset "Base.one(Milankovich)" begin
        # Test one() with default type
        mil_one = one(Milankovich)
        @test mil_one isa Milankovich{Float64}
        @test mil_one.hx == 0.0
        @test mil_one.hy == 0.0
        @test mil_one.hz == 0.0
        @test mil_one.ex == 0.0
        @test mil_one.ey == 0.0
        @test mil_one.ez == 0.0
        @test mil_one.L == 0.0
        
        # Test one() with specified type
        mil_one_f32 = one(Milankovich, T=Float32)
        @test mil_one_f32 isa Milankovich{Float32}
        @test mil_one_f32.hx == 0.0f0
    end
    
    @testset "Base.getindex - All Indices" begin
        # Test getindex for all 7 indices
        mil = Milankovich(1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 1.5)
        
        @test mil[1] == 1.0   # hx
        @test mil[2] == 2.0   # hy
        @test mil[3] == 3.0   # hz
        @test mil[4] == 0.1   # ex
        @test mil[5] == 0.2   # ey
        @test mil[6] == 0.3   # ez
        @test mil[7] == 1.5   # L
        
        # Test out of bounds
        @test_throws BoundsError mil[0]
        @test_throws BoundsError mil[8]
        @test_throws BoundsError mil[-1]
    end
    
    @testset "Transformations to/from Cartesian" begin
        μ = 398600.4418  # Earth gravitational parameter
        
        # Create a Milankovich coordinate
        mil = Milankovich(50000.0, 0.0, 50000.0, 0.05, 0.0, 0.0, 1.0)
        
        # Convert to Cartesian
        cart = Cartesian(mil, μ)
        @test cart isa Cartesian
        @test all(isfinite.(params(cart)))
        
        # Convert back to Milankovich
        mil_back = Milankovich(cart, μ)
        @test mil_back isa Milankovich
        @test all(isfinite.(params(mil_back)))
    end
    
    @testset "Round-Trip Transformations" begin
        μ = 398600.4418
        
        # Test round-trip: Milankovich → Cartesian → Milankovich
        mil_orig = Milankovich(50000.0, 0.0, 50000.0, 0.05, 0.0, 0.0, 1.0)
        cart = Cartesian(mil_orig, μ)
        mil_back = Milankovich(cart, μ)
        
        @test mil_back.hx ≈ mil_orig.hx atol=1e-6 rtol=1e-8
        @test mil_back.hy ≈ mil_orig.hy atol=1e-6 rtol=1e-8
        @test mil_back.hz ≈ mil_orig.hz atol=1e-6 rtol=1e-8
        @test mil_back.ex ≈ mil_orig.ex atol=1e-8 rtol=1e-8
        @test mil_back.ey ≈ mil_orig.ey atol=1e-8 rtol=1e-8
        @test mil_back.ez ≈ mil_orig.ez atol=1e-8 rtol=1e-8
    end
    
    @testset "Different Orbit Types" begin
        μ = 398600.4418
        
        # Test with circular orbit (small eccentricity)
        mil_circ = Milankovich(50000.0, 0.0, 50000.0, 0.001, 0.0, 0.0, 0.5)
        cart_circ = Cartesian(mil_circ, μ)
        @test all(isfinite.(params(cart_circ)))
        
        # Test with elliptical orbit
        mil_ellip = Milankovich(50000.0, 0.0, 45000.0, 0.2, 0.1, 0.0, 1.0)
        cart_ellip = Cartesian(mil_ellip, μ)
        @test all(isfinite.(params(cart_ellip)))
        
        # Test with inclined orbit (non-zero hz)
        mil_incl = Milankovich(50000.0, 5000.0, 48000.0, 0.1, 0.0, 0.05, 1.5)
        cart_incl = Cartesian(mil_incl, μ)
        @test all(isfinite.(params(cart_incl)))
    end
    
    @testset "Type Stability" begin
        μ = 398600.4418
        
        # Test Float64
        mil_f64 = Milankovich(50000.0, 0.0, 50000.0, 0.05, 0.0, 0.0, 1.0)
        cart_f64 = Cartesian(mil_f64, μ)
        @test eltype(cart_f64) === Float64
        
        # Test Float32
        mil_f32 = Milankovich(50000.0f0, 0.0f0, 50000.0f0, 0.05f0, 0.0f0, 0.0f0, 1.0f0)
        cart_f32 = Cartesian(mil_f32, Float32(μ))
        @test eltype(cart_f32) <: AbstractFloat
    end
    
    @testset "Edge Cases" begin
        μ = 398600.4418
        
        # Test with zero angular momentum components
        mil_zero_h = Milankovich(0.0, 0.0, 50000.0, 0.05, 0.0, 0.0, 1.0)
        cart_zero_h = Cartesian(mil_zero_h, μ)
        @test all(isfinite.(params(cart_zero_h)))
        
        # Test with zero eccentricity components
        mil_zero_e = Milankovich(50000.0, 0.0, 50000.0, 0.0, 0.0, 0.0, 1.0)
        cart_zero_e = Cartesian(mil_zero_e, μ)
        @test all(isfinite.(params(cart_zero_e)))
    end
end
