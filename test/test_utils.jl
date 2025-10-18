@testset "Utility Functions" begin
    @testset "angle_between_vectors - Parallel Vectors" begin
        # Two parallel vectors (angle = 0)
        v1 = SVector{3,Float64}(1.0, 0.0, 0.0)
        v2 = SVector{3,Float64}(2.0, 0.0, 0.0)  # Same direction, different magnitude

        angle = AstroCoords.angle_between_vectors(v1, v2)
        @test angle ≈ 0.0 atol=1e-10
    end

    @testset "angle_between_vectors - Perpendicular Vectors" begin
        # Two perpendicular vectors (angle = 90°)
        v1 = SVector{3,Float64}(1.0, 0.0, 0.0)
        v2 = SVector{3,Float64}(0.0, 1.0, 0.0)

        angle = AstroCoords.angle_between_vectors(v1, v2)
        @test angle ≈ π/2 atol=1e-10
    end

    @testset "angle_between_vectors - Opposite Vectors" begin
        # Two opposite vectors (angle = 180°)
        v1 = SVector{3,Float64}(1.0, 0.0, 0.0)
        v2 = SVector{3,Float64}(-1.0, 0.0, 0.0)

        angle = AstroCoords.angle_between_vectors(v1, v2)
        @test angle ≈ π atol=1e-10
    end

    @testset "angle_between_vectors - Arbitrary Angle" begin
        # Test 45° angle
        v1 = SVector{3,Float64}(1.0, 0.0, 0.0)
        v2 = normalize(SVector{3,Float64}(1.0, 1.0, 0.0))

        angle = AstroCoords.angle_between_vectors(v1, v2)
        @test angle ≈ π/4 atol=1e-10
    end

    @testset "angle_between_vectors - 3D Arbitrary Vectors" begin
        # Test with general 3D vectors
        v1 = SVector{3,Float64}(1.0, 2.0, 3.0)
        v2 = SVector{3,Float64}(4.0, 5.0, 6.0)

        angle = AstroCoords.angle_between_vectors(v1, v2)

        # Verify using dot product formula (as cross-check)
        cos_angle = dot(normalize(v1), normalize(v2))
        expected_angle = acos(clamp(cos_angle, -1.0, 1.0))
        @test angle ≈ expected_angle atol=1e-9
    end

    @testset "angle_between_vectors - Type Promotion (Float32/Float64)" begin
        # Mix Float32 and Float64 types
        v1 = SVector{3,Float32}(1.0f0, 0.0f0, 0.0f0)
        v2 = SVector{3,Float64}(0.0, 1.0, 0.0)

        angle = AstroCoords.angle_between_vectors(v1, v2)
        @test eltype(angle) == Float64  # Should promote to Float64
        @test angle ≈ π/2 atol=1e-6
    end

    @testset "angle_between_vectors - Non-Unit Vectors" begin
        # Test with non-normalized vectors (function should handle normalization)
        v1 = SVector{3,Float64}(5.0, 0.0, 0.0)
        v2 = SVector{3,Float64}(0.0, 3.0, 0.0)

        angle = AstroCoords.angle_between_vectors(v1, v2)
        @test angle ≈ π/2 atol=1e-10
    end

    @testset "angle_between_vectors - Different Magnitudes" begin
        # Vectors with very different magnitudes
        v1 = SVector{3,Float64}(1000.0, 0.0, 0.0)
        v2 = SVector{3,Float64}(0.001, 0.001, 0.0)

        angle = AstroCoords.angle_between_vectors(v1, v2)
        @test angle ≈ π/4 atol=1e-9
    end

    @testset "angle_between_vectors - Near-Zero Vector Stability" begin
        # Test numerical stability with small (but not zero) vectors
        v1 = SVector{3,Float64}(1e-10, 0.0, 0.0)
        v2 = SVector{3,Float64}(0.0, 1e-10, 0.0)

        angle = AstroCoords.angle_between_vectors(v1, v2)
        @test angle ≈ π/2 atol=1e-8  # Slightly relaxed tolerance for small vectors
    end

    @testset "angle_between_vectors - Acute vs Obtuse" begin
        # Test acute angle (< 90°)
        v1 = SVector{3,Float64}(1.0, 0.0, 0.0)
        v2 = SVector{3,Float64}(1.0, 0.5, 0.0)
        angle_acute = AstroCoords.angle_between_vectors(v1, v2)
        @test angle_acute < π/2
        @test angle_acute > 0.0

        # Test obtuse angle (> 90°)
        v3 = SVector{3,Float64}(-1.0, 0.5, 0.0)
        angle_obtuse = AstroCoords.angle_between_vectors(v1, v3)
        @test angle_obtuse > π/2
        @test angle_obtuse < π
    end

    @testset "angle_between_vectors - Multiple Random Cases" begin
        # Test multiple random vector pairs
        for _ in 1:20
            v1 = randn(SVector{3,Float64})
            v2 = randn(SVector{3,Float64})

            angle = AstroCoords.angle_between_vectors(v1, v2)

            # Angle should be in valid range [0, π]
            @test 0.0 <= angle <= π
            @test isfinite(angle)
        end
    end
end
