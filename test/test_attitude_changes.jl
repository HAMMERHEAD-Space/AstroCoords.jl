@testset "Attitude Transformations" begin
    @testset "EP2MRP - Identity Rotation" begin
        # Identity rotation: β = [1, 0, 0, 0]
        β = SVector{4,Float64}(1.0, 0.0, 0.0, 0.0)
        σ = AstroCoords.EP2MRP(β)

        # For identity, MRPs should be zero
        @test σ[1] ≈ 0.0 atol=1e-14
        @test σ[2] ≈ 0.0 atol=1e-14
        @test σ[3] ≈ 0.0 atol=1e-14
    end

    @testset "EP2MRP - 90° Rotation About X-axis" begin
        # 90° rotation about X-axis
        # β = [cos(θ/2), sin(θ/2)*axis]
        θ = π/2
        β = SVector{4,Float64}(cos(θ/2), sin(θ/2), 0.0, 0.0)
        σ = AstroCoords.EP2MRP(β)

        # σ = tan(θ/4) * axis
        expected_mag = tan(θ/4)
        @test σ[1] ≈ expected_mag atol=1e-14
        @test σ[2] ≈ 0.0 atol=1e-14
        @test σ[3] ≈ 0.0 atol=1e-14
    end

    @testset "EP2MRP - 90° Rotation About Y-axis" begin
        θ = π/2
        β = SVector{4,Float64}(cos(θ/2), 0.0, sin(θ/2), 0.0)
        σ = AstroCoords.EP2MRP(β)

        expected_mag = tan(θ/4)
        @test σ[1] ≈ 0.0 atol=1e-14
        @test σ[2] ≈ expected_mag atol=1e-14
        @test σ[3] ≈ 0.0 atol=1e-14
    end

    @testset "EP2MRP - 90° Rotation About Z-axis" begin
        θ = π/2
        β = SVector{4,Float64}(cos(θ/2), 0.0, 0.0, sin(θ/2))
        σ = AstroCoords.EP2MRP(β)

        expected_mag = tan(θ/4)
        @test σ[1] ≈ 0.0 atol=1e-14
        @test σ[2] ≈ 0.0 atol=1e-14
        @test σ[3] ≈ expected_mag atol=1e-14
    end

    @testset "EP2MRP - 180° Rotation (Edge Case)" begin
        # 180° rotation about X-axis: β = [0, 1, 0, 0]
        # This is a singularity case where β0 = 0
        β = SVector{4,Float64}(0.0, 1.0, 0.0, 0.0)
        σ = AstroCoords.EP2MRP(β)

        # At 180°, tan(θ/4) = tan(45°) = 1
        @test σ[1] ≈ 1.0 atol=1e-14
        @test σ[2] ≈ 0.0 atol=1e-14
        @test σ[3] ≈ 0.0 atol=1e-14
    end

    @testset "EP2MRP - Near-Singularity (β0 ≈ -1)" begin
        # Test numerical stability when β0 is close to -1
        # This would make denominator (1 + β0) very small
        β = SVector{4,Float64}(-0.999, 0.04, 0.02, 0.01)
        # Normalize to unit quaternion
        β = β / norm(β)
        σ = AstroCoords.EP2MRP(β)

        # Should still produce valid MRPs (check they're finite)
        @test isfinite(σ[1])
        @test isfinite(σ[2])
        @test isfinite(σ[3])
    end

    @testset "EP2MRP - Small Angle Rotation" begin
        # Small rotation (10 degrees) about arbitrary axis
        θ = deg2rad(10)
        axis = normalize(SVector{3,Float64}(1.0, 2.0, 3.0))
        β = SVector{4,Float64}(
            cos(θ/2), sin(θ/2)*axis[1], sin(θ/2)*axis[2], sin(θ/2)*axis[3]
        )
        σ = AstroCoords.EP2MRP(β)

        # For small angles, σ ≈ (θ/4) * axis
        expected_mag = tan(θ/4)
        @test norm(σ) ≈ expected_mag atol=1e-14
    end

    @testset "EP2MRP - Type Promotion" begin
        # Mix Float32 and Float64
        β = SVector{4}(Float32(0.7071), Float64(0.7071), Float32(0.0), Float64(0.0))
        σ = AstroCoords.EP2MRP(β)

        @test eltype(σ) == Float64  # Should promote to Float64
        @test isfinite(σ[1])
    end

    @testset "MRP2EP - Identity" begin
        # Zero MRP should give identity EP
        σ = SVector{3,Float64}(0.0, 0.0, 0.0)
        β = AstroCoords.MRP2EP(σ)

        @test β[1] ≈ 1.0 atol=1e-14 # β0 = 1 for identity
        @test β[2] ≈ 0.0 atol=1e-14
        @test β[3] ≈ 0.0 atol=1e-14
        @test β[4] ≈ 0.0 atol=1e-14
    end

    @testset "MRP2EP - 90° Rotation" begin
        # σ corresponding to 90° rotation about X
        θ = π/2
        σ = SVector{3,Float64}(tan(θ/4), 0.0, 0.0)
        β = AstroCoords.MRP2EP(σ)

        @test β[1] ≈ cos(θ/2) atol=1e-14
        @test β[2] ≈ sin(θ/2) atol=1e-14
        @test β[3] ≈ 0.0 atol=1e-14
        @test β[4] ≈ 0.0 atol=1e-14
    end

    @testset "MRP2EP - Round Trip" begin
        # Test EP → MRP → EP round trip
        θ = π/3  # 60 degrees
        axis = normalize(SVector{3,Float64}(1.0, 1.0, 1.0))
        β_orig = SVector{4,Float64}(
            cos(θ/2), sin(θ/2)*axis[1], sin(θ/2)*axis[2], sin(θ/2)*axis[3]
        )

        σ = AstroCoords.EP2MRP(β_orig)
        β_back = AstroCoords.MRP2EP(σ)

        # Should recover original EP (within tolerance)
        @test β_back[1] ≈ β_orig[1] atol=1e-14
        @test β_back[2] ≈ β_orig[2] atol=1e-14
        @test β_back[3] ≈ β_orig[3] atol=1e-14
        @test β_back[4] ≈ β_orig[4] atol=1e-14
    end

    @testset "MRP2EP - Multiple Round Trips" begin
        # Test multiple random rotations for bijection
        for _ in 1:10
            θ = π * rand()  # Random angle [0, π]
            axis = normalize(randn(SVector{3,Float64}))
            β_orig = SVector{4,Float64}(
                cos(θ/2), sin(θ/2)*axis[1], sin(θ/2)*axis[2], sin(θ/2)*axis[3]
            )

            σ = AstroCoords.EP2MRP(β_orig)
            β_back = AstroCoords.MRP2EP(σ)

            @test β_back[1] ≈ β_orig[1] atol=1e-14
            @test β_back[2] ≈ β_orig[2] atol=1e-14
            @test β_back[3] ≈ β_orig[3] atol=1e-14
            @test β_back[4] ≈ β_orig[4] atol=1e-14
        end
    end

    @testset "MRP2EP - Large MRP Values" begin
        # Test with larger MRP values (corresponding to larger rotations)
        σ = SVector{3,Float64}(2.0, 1.5, 1.0)
        β = AstroCoords.MRP2EP(σ)

        # EP should still be normalized
        @test norm(β) ≈ 1.0 atol=1e-14
    end
end
