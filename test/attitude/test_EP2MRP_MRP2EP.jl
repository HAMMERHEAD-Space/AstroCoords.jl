using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "EP2MRP and MRP2EP Conversions" begin
    @testset "EP2MRP: Standard Quaternion" begin
        # Tests lines 14-24
        # 90° rotation about x-axis: q = [cos(45°), sin(45°), 0, 0]
        β0 = cos(π/4)
        β1 = sin(π/4)
        β2 = 0.0
        β3 = 0.0
        β = SVector{4}(β0, β1, β2, β3)
        
        σ = EP2MRP(β)
        
        @test length(σ) == 3
        @test !any(isnan, σ)
        
        # For 90° rotation: MRP should be tan(45°) = 1.0 in the axis direction
        @test σ[1] ≈ tan(π/4) atol=1e-14
        @test σ[2] ≈ 0.0 atol=1e-14
        @test σ[3] ≈ 0.0 atol=1e-14
    end

    @testset "MRP2EP: Roundtrip Validation" begin
        # Tests lines 40-48
        β0_orig = 0.7071
        β1_orig = 0.7071
        β2_orig = 0.0
        β3_orig = 0.0
        β_orig = SVector{4}(β0_orig, β1_orig, β2_orig, β3_orig)
        
        σ = EP2MRP(β_orig)
        β_restored = MRP2EP(σ)
        
        # Quaternions can differ by sign, so check both
        @test (β_restored[1] ≈ β0_orig atol=1e-14 && 
               β_restored[2] ≈ β1_orig atol=1e-14 && 
               β_restored[3] ≈ β2_orig atol=1e-14 && 
               β_restored[4] ≈ β3_orig atol=1e-14) ||
              (β_restored[1] ≈ -β0_orig atol=1e-14 && 
               β_restored[2] ≈ -β1_orig atol=1e-14 && 
               β_restored[3] ≈ -β2_orig atol=1e-14 && 
               β_restored[4] ≈ -β3_orig atol=1e-14)
    end

    @testset "Edge Case: β0 ≈ -1 (Near-Singularity)" begin
        # Tests line 18 epsilon handling
        β0 = -1.0 + 1e-10  # Very close to -1
        β1 = 1e-5
        β2 = 1e-5
        β3 = 1e-5
        
        # Normalize to unit quaternion
        norm_factor = sqrt(β0^2 + β1^2 + β2^2 + β3^2)
        β = SVector{4}(β0, β1, β2, β3) ./ norm_factor
        
        σ = EP2MRP(β)
        
        # Should handle gracefully without division by zero
        @test !any(isnan, σ)
        @test !any(isinf, σ)
    end

    @testset "Multiple Rotation Examples" begin
        # 0° rotation (identity)
        β_0deg = SVector{4}(1.0, 0.0, 0.0, 0.0)
        σ_0deg = EP2MRP(β_0deg)
        @test norm(σ_0deg) ≈ 0.0 atol=1e-14
        
        # 180° rotation about x-axis: q = [0, 1, 0, 0]
        β_180x = SVector{4}(0.0, 1.0, 0.0, 0.0)
        σ_180x = EP2MRP(β_180x)
        # MRP becomes infinite at 180°, but should be very large
        @test abs(σ_180x[1]) > 100.0
        
        # 90° rotation about y-axis
        β_90y = SVector{4}(cos(π/4), 0.0, sin(π/4), 0.0)
        σ_90y = EP2MRP(β_90y)
        @test σ_90y[1] ≈ 0.0 atol=1e-14
        @test σ_90y[2] ≈ tan(π/4) atol=1e-14
        @test σ_90y[3] ≈ 0.0 atol=1e-14
        
        # 90° rotation about z-axis
        β_90z = SVector{4}(cos(π/4), 0.0, 0.0, sin(π/4))
        σ_90z = EP2MRP(β_90z)
        @test σ_90z[1] ≈ 0.0 atol=1e-14
        @test σ_90z[2] ≈ 0.0 atol=1e-14
        @test σ_90z[3] ≈ tan(π/4) atol=1e-14
    end

    @testset "Type Promotion" begin
        # Tests line 15 - type promotion
        β_f32 = SVector{4,Float32}(0.7071f0, 0.7071f0, 0.0f0, 0.0f0)
        σ_f32 = EP2MRP(β_f32)
        @test eltype(σ_f32) == Float32
        
        # MRP to EP promotion
        σ_f32_mrp = SVector{3,Float32}(1.0f0, 0.0f0, 0.0f0)
        β_f32_restored = MRP2EP(σ_f32_mrp)
        @test eltype(β_f32_restored) == Float32
    end

    @testset "Quaternion Normalization Check" begin
        # Ensure that after roundtrip, quaternion is normalized
        β_orig = SVector{4}(0.5, 0.5, 0.5, 0.5)
        σ = EP2MRP(β_orig)
        β_restored = MRP2EP(σ)
        
        # Check normalization
        norm_restored = sqrt(sum(β_restored .^ 2))
        @test norm_restored ≈ 1.0 atol=1e-14
    end

    @testset "Small Rotation Accuracy" begin
        # Small rotation (1 degree about x-axis)
        angle = π/180
        β_small = SVector{4}(cos(angle/2), sin(angle/2), 0.0, 0.0)
        σ_small = EP2MRP(β_small)
        β_small_restored = MRP2EP(σ_small)
        
        # Should preserve accuracy for small rotations
        @test (β_small_restored[1] ≈ β_small[1] atol=1e-14 &&
               β_small_restored[2] ≈ β_small[2] atol=1e-14) ||
              (β_small_restored[1] ≈ -β_small[1] atol=1e-14 &&
               β_small_restored[2] ≈ -β_small[2] atol=1e-14)
    end

    @testset "Arbitrary 3D Rotation" begin
        # Rotation about arbitrary axis [1,1,1]
        axis = normalize(SVector{3}(1.0, 1.0, 1.0))
        angle = π/3  # 60 degrees
        
        β0 = cos(angle/2)
        βv = sin(angle/2) * axis
        β = SVector{4}(β0, βv[1], βv[2], βv[3])
        
        σ = EP2MRP(β)
        β_restored = MRP2EP(σ)
        
        # Check roundtrip (accounting for sign ambiguity)
        norm_diff1 = norm(β - β_restored)
        norm_diff2 = norm(β + β_restored)
        @test min(norm_diff1, norm_diff2) < 1e-13
    end
end
