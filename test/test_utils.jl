using Test
using AstroCoords
using LinearAlgebra

@testset "Utility Functions" begin
    
    @testset "angle_between_vectors" begin
        #@lit {citation="Strang 2016", id="strang2016", ref="Linear Algebra, Section 1.2"}
        @testset "Parallel vectors (angle = 0)" begin
            v1 = [1.0, 0.0, 0.0]
            v2 = [2.0, 0.0, 0.0]  # Parallel, same direction
            
            angle = angle_between_vectors(v1, v2)
            
            @test angle ≈ 0.0 atol=1e-10
        end
        
        @testset "Anti-parallel vectors (angle = π)" begin
            v1 = [1.0, 0.0, 0.0]
            v2 = [-1.0, 0.0, 0.0]  # Anti-parallel
            
            angle = angle_between_vectors(v1, v2)
            
            @test angle ≈ π atol=1e-10
        end
        
        @testset "Orthogonal vectors (angle = π/2)" begin
            v1 = [1.0, 0.0, 0.0]
            v2 = [0.0, 1.0, 0.0]  # Orthogonal
            
            angle = angle_between_vectors(v1, v2)
            
            @test angle ≈ π/2 atol=1e-10
        end
        
        @testset "45° angle" begin
            v1 = [1.0, 0.0, 0.0]
            v2 = [1.0, 1.0, 0.0]  # 45° in xy-plane
            
            angle = angle_between_vectors(v1, v2)
            
            @test angle ≈ π/4 atol=1e-10
        end
        
        @testset "Arbitrary 3D vectors" begin
            v1 = [1.0, 2.0, 3.0]
            v2 = [4.0, 5.0, 6.0]
            
            angle = angle_between_vectors(v1, v2)
            
            # Compute expected angle manually
            cos_angle = dot(v1, v2) / (norm(v1) * norm(v2))
            expected = acos(clamp(cos_angle, -1.0, 1.0))
            
            @test angle ≈ expected atol=1e-10 rtol=1e-10
            @test 0.0 <= angle <= π
        end
        
        @testset "Numerical stability: nearly parallel" begin
            # Vectors that are nearly parallel (dot ≈ 1)
            v1 = [1.0, 0.0, 0.0]
            v2 = [1.0, 1e-10, 0.0]
            
            angle = angle_between_vectors(v1, v2)
            
            # Should be very small but not exactly zero
            @test angle ≈ 0.0 atol=1e-8
            @test angle >= 0.0
        end
        
        @testset "Numerical stability: nearly anti-parallel" begin
            # Vectors that are nearly anti-parallel (dot ≈ -1)
            v1 = [1.0, 0.0, 0.0]
            v2 = [-1.0, 1e-10, 0.0]
            
            angle = angle_between_vectors(v1, v2)
            
            # Should be close to π
            @test angle ≈ π atol=1e-8
            @test angle <= π
        end
        
        @testset "Zero-length vector (should error)" begin
            v1 = [0.0, 0.0, 0.0]
            v2 = [1.0, 0.0, 0.0]
            
            # normalize() should throw an error for zero vector
            @test_throws Exception angle_between_vectors(v1, v2)
        end
        
        @testset "Mixed precision: Float32 and Float64" begin
            v1_f32 = Float32[1.0, 2.0, 3.0]
            v2_f64 = Float64[4.0, 5.0, 6.0]
            
            # Should use promote_type to handle mixed types
            angle = angle_between_vectors(v1_f32, v2_f64)
            
            # Result should be Float64 (promoted type)
            @test typeof(angle) == Float64
            @test isfinite(angle)
            @test 0.0 <= angle <= π
        end
        
        @testset "Different vector lengths" begin
            v1 = [1.0, 0.0, 0.0]
            v2 = [100.0, 0.0, 0.0]
            
            angle = angle_between_vectors(v1, v2)
            
            # Angle should be independent of magnitude
            @test angle ≈ 0.0 atol=1e-10
        end
        
        @testset "Symmetry property" begin
            v1 = [1.0, 2.0, 3.0]
            v2 = [4.0, 5.0, 6.0]
            
            angle1 = angle_between_vectors(v1, v2)
            angle2 = angle_between_vectors(v2, v1)
            
            # Should be symmetric
            @test angle1 ≈ angle2 atol=1e-10
        end
        
        @testset "Clamping behavior" begin
            # Test that clamping prevents acos domain errors
            # This is implicitly tested by the "nearly parallel" cases,
            # but we can also test the dot product directly
            
            v1 = normalize([1.0, 0.0, 0.0])
            v2 = normalize([1.0 + 1e-15, 0.0, 0.0])  # Numerical error might make dot > 1
            
            angle = angle_between_vectors(v1, v2)
            
            # Should not throw domain error from acos
            @test isfinite(angle)
            @test angle >= 0.0
            @test angle <= π
        end
    end
end
