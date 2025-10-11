using Test
using AstroCoords
using LinearAlgebra

@testset "Utility Functions" begin
    @testset "angle_between_vectors" begin
        @testset "Orthogonal vectors (90°)" begin
            v1 = [1.0, 0.0, 0.0]
            v2 = [0.0, 1.0, 0.0]
            
            angle = angle_between_vectors(v1, v2)
            @test angle ≈ π/2 atol=1e-15
        end
        
        @testset "Parallel vectors (0°)" begin
            v1 = [1.0, 2.0, 3.0]
            v2 = [2.0, 4.0, 6.0]  # Parallel (scaled)
            
            angle = angle_between_vectors(v1, v2)
            @test angle ≈ 0.0 atol=1e-15
        end
        
        @testset "Antiparallel vectors (180°)" begin
            v1 = [1.0, 2.0, 3.0]
            v2 = [-1.0, -2.0, -3.0]  # Antiparallel
            
            angle = angle_between_vectors(v1, v2)
            @test angle ≈ π atol=1e-15
        end
        
        @testset "Known angles: 30°, 45°, 60°" begin
            # 30° angle between [1,0,0] and [cos(30°), sin(30°), 0]
            v1 = [1.0, 0.0, 0.0]
            v2 = [cos(π/6), sin(π/6), 0.0]
            angle = angle_between_vectors(v1, v2)
            @test angle ≈ π/6 atol=1e-15
            
            # 45° angle
            v3 = [cos(π/4), sin(π/4), 0.0]
            angle = angle_between_vectors(v1, v3)
            @test angle ≈ π/4 atol=1e-15
            
            # 60° angle
            v4 = [cos(π/3), sin(π/3), 0.0]
            angle = angle_between_vectors(v1, v4)
            @test angle ≈ π/3 atol=1e-15
        end
        
        @testset "3D vectors at known angles" begin
            # Angle between [1,1,0] and [1,0,1]
            v1 = [1.0, 1.0, 0.0]
            v2 = [1.0, 0.0, 1.0]
            
            # cos(θ) = (v1·v2)/(|v1||v2|) = 1/(√2·√2) = 1/2 → θ = 60°
            angle = angle_between_vectors(v1, v2)
            @test angle ≈ π/3 atol=1e-15
        end
        
        @testset "Unnormalized vectors" begin
            # Should work with unnormalized vectors of different magnitudes
            v1 = [100.0, 0.0, 0.0]
            v2 = [0.0, 0.001, 0.0]
            
            angle = angle_between_vectors(v1, v2)
            @test angle ≈ π/2 atol=1e-15
        end
        
        @testset "Numerical stability: nearly parallel" begin
            # Vectors with very small angle (1e-12 radians ≈ 5.7e-11 degrees)
            θ = 1e-12
            v1 = [1.0, 0.0, 0.0]
            v2 = [cos(θ), sin(θ), 0.0]
            
            angle = angle_between_vectors(v1, v2)
            # Should be accurate even for very small angles
            @test angle ≈ θ atol=1e-14 rtol=1e-10
        end
        
        @testset "Numerical stability: nearly antiparallel" begin
            # Vectors with angle very close to π
            θ = π - 1e-12
            v1 = [1.0, 0.0, 0.0]
            v2 = [cos(θ), sin(θ), 0.0]
            
            angle = angle_between_vectors(v1, v2)
            # Should be accurate even near π
            @test angle ≈ θ atol=1e-14 rtol=1e-10
        end
        
        @testset "Zero-length vector handling" begin
            v1 = [1.0, 0.0, 0.0]
            v2 = [0.0, 0.0, 0.0]  # Zero vector
            
            # Should either throw an error or return NaN/Inf
            # Testing that it doesn't silently give wrong answer
            result = angle_between_vectors(v1, v2)
            @test isnan(result) || isinf(result) || result isa Exception
        end
        
        @testset "Symmetry property" begin
            # angle(v1, v2) should equal angle(v2, v1)
            v1 = [1.0, 2.0, 3.0]
            v2 = [4.0, -1.0, 2.0]
            
            angle1 = angle_between_vectors(v1, v2)
            angle2 = angle_between_vectors(v2, v1)
            
            @test angle1 ≈ angle2 atol=1e-15
        end
        
        @testset "Range: angles should be in [0, π]" begin
            # All angles between vectors should be in [0, π]
            test_vectors = [
                ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
                ([1.0, 1.0, 1.0], [1.0, -1.0, 0.0]),
                ([-1.0, 2.0, -3.0], [0.5, 0.5, 0.5]),
            ]
            
            for (v1, v2) in test_vectors
                angle = angle_between_vectors(v1, v2)
                @test 0.0 <= angle <= π
            end
        end
    end
end
