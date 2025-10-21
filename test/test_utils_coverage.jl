using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "Utility Functions Coverage" begin
    @testset "angle_between_vectors - parallel vectors" begin
        # Same direction
        v1 = [1.0, 0.0, 0.0]
        v2 = [2.0, 0.0, 0.0]
        angle = angle_between_vectors(v1, v2)
        @test angle ≈ 0.0 atol=1e-14
        
        # Same direction, different magnitudes
        v3 = [0.0, 5.0, 0.0]
        v4 = [0.0, 0.1, 0.0]
        angle2 = angle_between_vectors(v3, v4)
        @test angle2 ≈ 0.0 atol=1e-14
        
        # SVector inputs
        sv1 = SVector{3}(1.0, 2.0, 3.0)
        sv2 = SVector{3}(2.0, 4.0, 6.0)
        angle_sv = angle_between_vectors(sv1, sv2)
        @test angle_sv ≈ 0.0 atol=1e-14
    end
    
    @testset "angle_between_vectors - antiparallel vectors" begin
        # Opposite directions
        v1 = [1.0, 0.0, 0.0]
        v2 = [-1.0, 0.0, 0.0]
        angle = angle_between_vectors(v1, v2)
        @test angle ≈ π atol=1e-14
        
        # Different magnitudes
        v3 = [0.0, 3.0, 0.0]
        v4 = [0.0, -7.0, 0.0]
        angle2 = angle_between_vectors(v3, v4)
        @test angle2 ≈ π atol=1e-14
        
        # 3D case
        v5 = [1.0, 2.0, 3.0]
        v6 = [-2.0, -4.0, -6.0]
        angle3 = angle_between_vectors(v5, v6)
        @test angle3 ≈ π atol=1e-14
    end
    
    @testset "angle_between_vectors - perpendicular vectors" begin
        # X and Y axes
        v1 = [1.0, 0.0, 0.0]
        v2 = [0.0, 1.0, 0.0]
        angle = angle_between_vectors(v1, v2)
        @test angle ≈ π/2 atol=1e-14
        
        # Y and Z axes
        v3 = [0.0, 2.0, 0.0]
        v4 = [0.0, 0.0, 3.0]
        angle2 = angle_between_vectors(v3, v4)
        @test angle2 ≈ π/2 atol=1e-14
        
        # General perpendicular
        v5 = [1.0, 1.0, 0.0]
        v6 = [-1.0, 1.0, 0.0]
        angle3 = angle_between_vectors(v5, v6)
        @test angle3 ≈ π/2 atol=1e-14
    end
    
    @testset "angle_between_vectors - arbitrary angles" begin
        # 45 degrees
        v1 = [1.0, 0.0, 0.0]
        v2 = [1.0, 1.0, 0.0]
        angle = angle_between_vectors(v1, v2)
        @test angle ≈ π/4 atol=1e-14
        
        # 60 degrees
        v3 = [1.0, 0.0, 0.0]
        v4 = [0.5, sqrt(3)/2, 0.0]
        angle2 = angle_between_vectors(v3, v4)
        @test angle2 ≈ π/3 atol=1e-14
        
        # 135 degrees
        v5 = [1.0, 0.0, 0.0]
        v6 = [-1.0, 1.0, 0.0]
        angle3 = angle_between_vectors(v5, v6)
        @test angle3 ≈ 3π/4 atol=1e-14
    end
    
    @testset "angle_between_vectors - numerical stability near 0°" begin
        # Nearly parallel vectors (small angle)
        v1 = [1.0, 0.0, 0.0]
        v2 = [1.0, 1e-10, 0.0]
        angle = angle_between_vectors(v1, v2)
        @test angle ≈ 1e-10 atol=1e-14 rtol=1e-3
        @test angle >= 0.0
        
        # Very close vectors
        v3 = [1.0, 1.0, 1.0]
        v4 = [1.0 + 1e-12, 1.0, 1.0]
        angle2 = angle_between_vectors(v3, v4)
        @test angle2 >= 0.0
        @test angle2 < 1e-10
    end
    
    @testset "angle_between_vectors - numerical stability near 180°" begin
        # Nearly antiparallel vectors
        v1 = [1.0, 0.0, 0.0]
        v2 = [-1.0, 1e-10, 0.0]
        angle = angle_between_vectors(v1, v2)
        @test angle ≈ π atol=1e-10
        @test angle <= π
        
        # Very close to 180°
        v3 = [1.0, 1.0, 1.0]
        v4 = [-1.0 - 1e-12, -1.0, -1.0]
        angle2 = angle_between_vectors(v3, v4)
        @test angle2 ≈ π atol=1e-10
        @test angle2 <= π
    end
    
    @testset "angle_between_vectors - type promotion" begin
        # Float32 and Float64
        v1_f32 = Float32[1.0, 0.0, 0.0]
        v2_f64 = Float64[0.0, 1.0, 0.0]
        angle = angle_between_vectors(v1_f32, v2_f64)
        @test angle isa Float64
        @test angle ≈ π/2 atol=1e-6
        
        # Int and Float64
        v1_int = [1, 0, 0]
        v2_float = [0.0, 1.0, 0.0]
        angle2 = angle_between_vectors(v1_int, v2_float)
        @test angle2 isa Float64
        @test angle2 ≈ π/2 atol=1e-14
    end
    
    @testset "angle_between_vectors - different vector types" begin
        # Regular Array
        v1_arr = [1.0, 0.0, 0.0]
        v2_arr = [0.0, 1.0, 0.0]
        angle_arr = angle_between_vectors(v1_arr, v2_arr)
        @test angle_arr ≈ π/2 atol=1e-14
        
        # SVector
        v1_sv = SVector{3}(1.0, 0.0, 0.0)
        v2_sv = SVector{3}(0.0, 1.0, 0.0)
        angle_sv = angle_between_vectors(v1_sv, v2_sv)
        @test angle_sv ≈ π/2 atol=1e-14
        
        # Mixed types
        angle_mixed = angle_between_vectors(v1_arr, v2_sv)
        @test angle_mixed ≈ π/2 atol=1e-14
    end
    
    @testset "angle_between_vectors - edge cases with very small magnitudes" begin
        # Small magnitude vectors (but not zero)
        v1 = [1e-8, 0.0, 0.0]
        v2 = [0.0, 1e-8, 0.0]
        angle = angle_between_vectors(v1, v2)
        @test angle ≈ π/2 atol=1e-10
        
        # One small, one normal
        v3 = [1e-10, 0.0, 0.0]
        v4 = [0.0, 1.0, 0.0]
        angle2 = angle_between_vectors(v3, v4)
        @test angle2 ≈ π/2 atol=1e-10
    end
    
    @testset "angle_between_vectors - 3D space rotation consistency" begin
        # Vectors in XY plane
        v1 = [1.0, 0.0, 0.0]
        v2 = [cos(π/6), sin(π/6), 0.0]
        angle = angle_between_vectors(v1, v2)
        @test angle ≈ π/6 atol=1e-14
        
        # Vectors in XZ plane
        v3 = [1.0, 0.0, 0.0]
        v4 = [cos(π/3), 0.0, sin(π/3)]
        angle2 = angle_between_vectors(v3, v4)
        @test angle2 ≈ π/3 atol=1e-14
        
        # General 3D
        v5 = [1.0, 1.0, 1.0]
        v6 = [1.0, 0.0, 0.0]
        angle3 = angle_between_vectors(v5, v6)
        expected = acos(1.0 / sqrt(3))
        @test angle3 ≈ expected atol=1e-14
    end
end
