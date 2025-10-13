using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "angle_between_vectors" begin
    @testset "Orthogonal Vectors (90°)" begin
        # Tests lines 16-22
        v1 = SVector{3}(1.0, 0.0, 0.0)
        v2 = SVector{3}(0.0, 1.0, 0.0)
        
        angle = angle_between_vectors(v1, v2)
        @test angle ≈ π/2 atol=1e-14
    end

    @testset "Parallel Vectors (0°)" begin
        # Tests lines 24-27
        v1 = SVector{3}(1.0, 0.0, 0.0)
        v2 = SVector{3}(2.0, 0.0, 0.0)
        
        angle = angle_between_vectors(v1, v2)
        @test angle ≈ 0.0 atol=1e-14
    end

    @testset "Antiparallel Vectors (180°)" begin
        # Tests lines 29-32
        v1 = SVector{3}(1.0, 0.0, 0.0)
        v2 = SVector{3}(-1.0, 0.0, 0.0)
        
        angle = angle_between_vectors(v1, v2)
        @test angle ≈ π atol=1e-14
    end

    @testset "Arbitrary Angles" begin
        # Tests lines 19-35 for general case
        
        # 30 degrees
        v1 = SVector{3}(1.0, 0.0, 0.0)
        v2 = SVector{3}(cos(π/6), sin(π/6), 0.0)
        angle_30 = angle_between_vectors(v1, v2)
        @test angle_30 ≈ π/6 atol=1e-14
        
        # 60 degrees
        v2_60 = SVector{3}(cos(π/3), sin(π/3), 0.0)
        angle_60 = angle_between_vectors(v1, v2_60)
        @test angle_60 ≈ π/3 atol=1e-14
        
        # 120 degrees
        v2_120 = SVector{3}(cos(2*π/3), sin(2*π/3), 0.0)
        angle_120 = angle_between_vectors(v1, v2_120)
        @test angle_120 ≈ 2*π/3 atol=1e-14
    end

    @testset "Type Promotion" begin
        # Tests line 19 - Vector{Int} vs Vector{Float64}
        v1_int = SVector{3,Int}(1, 0, 0)
        v2_float = SVector{3,Float64}(0.0, 1.0, 0.0)
        
        angle = angle_between_vectors(v1_int, v2_float)
        @test angle isa Float64
        @test angle ≈ π/2 atol=1e-14
        
        # Both Int
        v1_int2 = SVector{3,Int}(1, 0, 0)
        v2_int2 = SVector{3,Int}(0, 1, 0)
        angle_int = angle_between_vectors(v1_int2, v2_int2)
        @test angle_int isa Float64  # Should promote to Float64
        @test angle_int ≈ π/2 atol=1e-14
    end

    @testset "Numerical Stability" begin
        # Test the numerically stable formulation
        # Compare with dot product method (less stable for small/large angles)
        
        # Very small angle
        v1 = SVector{3}(1.0, 0.0, 0.0)
        v2 = SVector{3}(cos(1e-10), sin(1e-10), 0.0)
        angle_small = angle_between_vectors(v1, v2)
        @test angle_small ≈ 1e-10 rtol=1e-8  # Should handle small angles accurately
        
        # Angle close to 180°
        v2_near180 = SVector{3}(-cos(1e-10), -sin(1e-10), 0.0)
        angle_near180 = angle_between_vectors(v1, v2_near180)
        @test angle_near180 ≈ π rtol=1e-8  # Should handle angles near π accurately
    end

    @testset "Non-Unit Vectors" begin
        # Function should normalize internally (line 21-22)
        v1 = SVector{3}(5.0, 0.0, 0.0)
        v2 = SVector{3}(0.0, 3.0, 0.0)
        
        angle = angle_between_vectors(v1, v2)
        @test angle ≈ π/2 atol=1e-14  # Should still be 90° after normalization
    end

    @testset "3D Vectors" begin
        # Test with vectors not in xy-plane
        v1 = SVector{3}(1.0, 1.0, 0.0)
        v2 = SVector{3}(0.0, 1.0, 1.0)
        
        # Expected angle: arccos(1/sqrt(2) * 1/sqrt(2)) = arccos(1/2) = π/3
        angle = angle_between_vectors(v1, v2)
        expected = acos(dot(normalize(v1), normalize(v2)))
        @test angle ≈ expected atol=1e-14
    end

    @testset "Edge Case: Near-Zero Vectors" begin
        # Very small magnitude vectors
        v1 = SVector{3}(1e-20, 0.0, 0.0)
        v2 = SVector{3}(0.0, 1e-20, 0.0)
        
        # Should still compute angle correctly after normalization
        angle = angle_between_vectors(v1, v2)
        @test angle ≈ π/2 atol=1e-10
    end

    @testset "Angle Clamping Behavior" begin
        # Tests lines 30-35 - clamping to [0, π]
        v1 = SVector{3}(1.0, 0.0, 0.0)
        v2 = SVector{3}(1.0, 0.0, 0.0)
        
        angle = angle_between_vectors(v1, v2)
        @test 0.0 <= angle <= π
        @test angle ≈ 0.0 atol=1e-14
    end
end
