using Test
using AstroCoords
using StaticArrays

@testset "Milankovich Coverage Tests" begin
    @testset "Constructor: AbstractVector{T}" begin
        # Test constructor from generic AbstractVector
        vec = [1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 1.5]
        mil = Milankovich(vec)
        @test mil.hx ≈ 1.0 atol=1e-15
        @test mil.hy ≈ 2.0 atol=1e-15
        @test mil.hz ≈ 3.0 atol=1e-15
        @test mil.ex ≈ 0.1 atol=1e-15
        @test mil.ey ≈ 0.2 atol=1e-15
        @test mil.ez ≈ 0.3 atol=1e-15
        @test mil.L ≈ 1.5 atol=1e-15
        @test typeof(mil) === Milankovich{Float64}
    end

    @testset "Constructor: Individual args" begin
        # Test constructor with individual arguments
        mil = Milankovich(1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 1.5)
        @test mil.hx ≈ 1.0 atol=1e-15
        @test mil.hy ≈ 2.0 atol=1e-15
        @test mil.hz ≈ 3.0 atol=1e-15
        @test mil.ex ≈ 0.1 atol=1e-15
        @test mil.ey ≈ 0.2 atol=1e-15
        @test mil.ez ≈ 0.3 atol=1e-15
        @test mil.L ≈ 1.5 atol=1e-15
        
        # Test type promotion with mixed types
        mil_mixed = Milankovich(1, 2.0f0, 3.0, 0.1, 0.2, 0.3, 1.5)
        @test typeof(mil_mixed) === Milankovich{Float64}
    end

    @testset "Constructor: StaticVector" begin
        # Test constructor from StaticVector (SVector)
        svec = SVector{7}(1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 1.5)
        mil = Milankovich(svec)
        @test mil.hx ≈ 1.0 atol=1e-15
        @test mil.hy ≈ 2.0 atol=1e-15
        @test mil.hz ≈ 3.0 atol=1e-15
        @test mil.ex ≈ 0.1 atol=1e-15
        @test mil.ey ≈ 0.2 atol=1e-15
        @test mil.ez ≈ 0.3 atol=1e-15
        @test mil.L ≈ 1.5 atol=1e-15
        
        # Test Milankovich{T}(::StaticVector) constructor
        mil_typed = Milankovich{Float64}(svec)
        @test mil_typed.hx ≈ 1.0 atol=1e-15
        @test typeof(mil_typed) === Milankovich{Float64}
    end

    @testset "Edge case: Zero eccentricity (ex=ey=0, circular orbit)" begin
        # Circular orbit case: ex=ey=0, ez can be non-zero
        mil = Milankovich(1.5, 2.0, 3.5, 0.0, 0.0, 0.1, 2.0)
        @test mil.ex ≈ 0.0 atol=1e-15
        @test mil.ey ≈ 0.0 atol=1e-15
        @test mil.ez ≈ 0.1 atol=1e-15
        
        # Verify params() extraction works
        p = params(mil)
        @test p[4] ≈ 0.0 atol=1e-15  # ex
        @test p[5] ≈ 0.0 atol=1e-15  # ey
    end

    @testset "Edge case: Equatorial orbit (hx=hy=0, i≈0)" begin
        # Equatorial orbit: hx=hy=0, hz non-zero
        mil = Milankovich(0.0, 0.0, 4.0, 0.2, 0.1, 0.15, 1.0)
        @test mil.hx ≈ 0.0 atol=1e-15
        @test mil.hy ≈ 0.0 atol=1e-15
        @test mil.hz ≈ 4.0 atol=1e-15
        
        # Verify params() extraction
        p = params(mil)
        @test p[1] ≈ 0.0 atol=1e-15  # hx
        @test p[2] ≈ 0.0 atol=1e-15  # hy
        @test p[3] ≈ 4.0 atol=1e-15  # hz
    end

    @testset "params() function" begin
        # Test params() returns correct SVector
        mil = Milankovich(1.1, 2.2, 3.3, 0.11, 0.22, 0.33, 1.234)
        p = params(mil)
        @test typeof(p) === SVector{7, Float64}
        @test p[1] ≈ 1.1 atol=1e-15
        @test p[2] ≈ 2.2 atol=1e-15
        @test p[3] ≈ 3.3 atol=1e-15
        @test p[4] ≈ 0.11 atol=1e-15
        @test p[5] ≈ 0.22 atol=1e-15
        @test p[6] ≈ 0.33 atol=1e-15
        @test p[7] ≈ 1.234 atol=1e-15
    end

    @testset "Type conversion edge cases" begin
        # Test Float32 input
        mil_f32 = Milankovich(Float32[1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 1.5])
        @test typeof(mil_f32) === Milankovich{Float32}
        
        # Test mixed precision promotion
        mil_mixed = Milankovich(1.0f0, 2.0, 3, 0.1f0, 0.2, 0.3f0, 1.5)
        @test typeof(mil_mixed) === Milankovich{Float64}
    end
end
