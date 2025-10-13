using Test
using AstroCoords
using LinearAlgebra

@testset "Utility Functions" begin
    @testset "angle_between_vectors - Orthogonal" begin
        # Orthogonal vectors should have angle π/2
        v1 = [1.0, 0.0, 0.0]
        v2 = [0.0, 1.0, 0.0]
        
        angle = AstroCoords.angle_between_vectors(v1, v2)
        
        @test angle ≈ π/2 atol=1e-14
    end

    @testset "angle_between_vectors - Parallel" begin
        # Parallel vectors should have angle 0
        v1 = [1.0, 2.0, 3.0]
        v2 = [2.0, 4.0, 6.0]
        
        angle = AstroCoords.angle_between_vectors(v1, v2)
        
        @test angle ≈ 0.0 atol=1e-14
    end

    @testset "angle_between_vectors - Antiparallel" begin
        # Antiparallel vectors should have angle π
        v1 = [1.0, 0.0, 0.0]
        v2 = [-1.0, 0.0, 0.0]
        
        angle = AstroCoords.angle_between_vectors(v1, v2)
        
        @test angle ≈ π atol=1e-14
    end

    @testset "angle_between_vectors - Arbitrary Non-Unit" begin
        # Non-unit vectors should be normalized internally
        v1 = [3.0, 4.0, 0.0]  # magnitude 5
        v2 = [0.0, 5.0, 12.0]  # magnitude 13
        
        # Manual calculation: after normalization
        # v1_norm = [0.6, 0.8, 0]
        # v2_norm = [0, 5/13, 12/13]
        # cos(θ) = v1_norm · v2_norm = 0.8 * (5/13) = 4/13
        expected_angle = acos(4/13)
        
        angle = AstroCoords.angle_between_vectors(v1, v2)
        
        @test angle ≈ expected_angle atol=1e-14
    end

    @testset "angle_between_vectors - 2D Vectors" begin
        # Test with 2D vectors (should still work)
        v1 = [1.0, 0.0]
        v2 = [0.0, 1.0]
        
        angle = AstroCoords.angle_between_vectors(v1, v2)
        
        @test angle ≈ π/2 atol=1e-14
    end

    @testset "angle_between_vectors - Known Angle (60°)" begin
        # Construct vectors with known 60° angle
        # v1 = [1, 0, 0], v2 = [cos(60°), sin(60°), 0] = [0.5, √3/2, 0]
        v1 = [1.0, 0.0, 0.0]
        v2 = [0.5, sqrt(3)/2, 0.0]
        
        angle = AstroCoords.angle_between_vectors(v1, v2)
        
        @test angle ≈ π/3 atol=1e-14
    end

    @testset "angle_between_vectors - Nearly Parallel" begin
        # Test numerical precision with nearly parallel vectors
        # cos(θ) ≈ 1 - ε for small ε
        v1 = [1.0, 0.0, 0.0]
        v2 = [1.0, 1e-10, 0.0]
        
        angle = AstroCoords.angle_between_vectors(v1, v2)
        
        # Should be a very small angle
        @test angle < 1e-9
        @test angle >= 0.0
    end

    @testset "angle_between_vectors - Type Promotion" begin
        # Test type promotion Float32 + Float64
        v1 = Float32[1.0, 0.0, 0.0]
        v2 = Float64[0.0, 1.0, 0.0]
        
        angle = AstroCoords.angle_between_vectors(v1, v2)
        
        @test angle ≈ π/2 atol=1e-6  # Float32 precision
        @test typeof(angle) == Float64
    end

    @testset "angle_between_vectors - Same Vector" begin
        # Same vector should give angle 0
        v1 = [1.0, 2.0, 3.0]
        
        angle = AstroCoords.angle_between_vectors(v1, v1)
        
        @test angle ≈ 0.0 atol=1e-14
    end

    @testset "angle_between_vectors - Coordinate Geometry" begin
        # Test with standard basis vectors
        ex = [1.0, 0.0, 0.0]
        ey = [0.0, 1.0, 0.0]
        ez = [0.0, 0.0, 1.0]
        
        @test AstroCoords.angle_between_vectors(ex, ey) ≈ π/2 atol=1e-14
        @test AstroCoords.angle_between_vectors(ey, ez) ≈ π/2 atol=1e-14
        @test AstroCoords.angle_between_vectors(ez, ex) ≈ π/2 atol=1e-14
    end
end
