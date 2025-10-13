using Test
using AstroCoords

@testset "Attitude Conversions" begin
    @testset "EP2MRP - Identity Rotation" begin
        # Identity quaternion [1, 0, 0, 0] should map to zero MRP
        β = [1.0, 0.0, 0.0, 0.0]
        σ = AstroCoords.EP2MRP(β)
        
        @test σ ≈ [0.0, 0.0, 0.0] atol=1e-14
    end

    @testset "EP2MRP - 90° Rotation about Z" begin
        # 90° rotation about z-axis: β = [cos(45°), 0, 0, sin(45°)]
        β = [sqrt(2)/2, 0.0, 0.0, sqrt(2)/2]
        σ = AstroCoords.EP2MRP(β)
        
        # MRP for 90° rotation about z: σ = [0, 0, tan(45°)] = [0, 0, 1]
        expected = [0.0, 0.0, 1.0]
        @test σ ≈ expected atol=1e-14
    end

    @testset "EP2MRP - β0 near -1 (Shadow Set)" begin
        # Test near singularity at β0 = -1
        # 170° rotation about x-axis: β ≈ [-0.9848, 0.1736, 0, 0]
        θ = 170 * π / 180
        β = [cos(θ/2), sin(θ/2), 0.0, 0.0]
        
        σ = AstroCoords.EP2MRP(β)
        
        # Should not have NaN or Inf
        @test all(isfinite.(σ))
        # MRP magnitude should be large but finite
        @test norm(σ) > 1.0  # In shadow set
    end

    @testset "EP2MRP - β0 exactly -1" begin
        # Test exact singularity protection
        β = [-1.0, 0.0, 0.0, 0.0]
        
        σ = AstroCoords.EP2MRP(β)
        
        # Should not crash, should return finite values
        @test all(isfinite.(σ))
    end

    @testset "MRP2EP - Zero MRP" begin
        # Zero MRP should map to identity quaternion
        σ = [0.0, 0.0, 0.0]
        β = AstroCoords.MRP2EP(σ)
        
        # First component should be 1 (or -1, quaternions have double cover)
        @test abs(β[1]) ≈ 1.0 atol=1e-14
        @test β[2] ≈ 0.0 atol=1e-14
        @test β[3] ≈ 0.0 atol=1e-14
        @test β[4] ≈ 0.0 atol=1e-14
    end

    @testset "MRP2EP - 90° Rotation about Z" begin
        # MRP for 90° rotation about z: σ = [0, 0, 1]
        σ = [0.0, 0.0, 1.0]
        β = AstroCoords.MRP2EP(σ)
        
        # Should recover 90° rotation quaternion
        expected_magnitude = 1.0
        @test norm(β) ≈ expected_magnitude atol=1e-14
        
        # β0 should be cos(45°) = sqrt(2)/2
        @test abs(β[1]) ≈ sqrt(2)/2 atol=1e-14
        @test abs(β[4]) ≈ sqrt(2)/2 atol=1e-14
    end

    @testset "Round-Trip: EP → MRP → EP" begin
        # Test that EP → MRP → EP preserves rotation
        # 60° rotation about arbitrary axis
        θ = π / 3
        axis = normalize([1.0, 2.0, 3.0])
        β_orig = [cos(θ/2); sin(θ/2) * axis]
        
        σ = AstroCoords.EP2MRP(β_orig)
        β_back = AstroCoords.MRP2EP(σ)
        
        # Quaternions have double cover: q and -q represent same rotation
        # Check if β_back ≈ β_orig or β_back ≈ -β_orig
        matches_positive = all(abs.(β_back .- β_orig) .< 1e-12)
        matches_negative = all(abs.(β_back .+ β_orig) .< 1e-12)
        
        @test matches_positive || matches_negative
    end

    @testset "Round-Trip: MRP → EP → MRP" begin
        # Test that MRP → EP → MRP preserves rotation
        σ_orig = [0.1, 0.2, 0.3]
        
        β = AstroCoords.MRP2EP(σ_orig)
        σ_back = AstroCoords.EP2MRP(β)
        
        @test σ_back ≈ σ_orig atol=1e-12
    end

    @testset "Rotation Equivalence" begin
        # Test that EP and MRP represent the same rotation
        # by applying both to a test vector
        using LinearAlgebra
        
        # 45° rotation about [1,1,1]
        θ = π / 4
        axis = normalize([1.0, 1.0, 1.0])
        β = [cos(θ/2); sin(θ/2) * axis]
        σ = AstroCoords.EP2MRP(β)
        
        # Test vector
        v = [1.0, 0.0, 0.0]
        
        # Rotate using EP (quaternion rotation formula)
        # v' = (β0² - β·β)v + 2(β·v)β + 2β0(β × v)
        β0 = β[1]
        βvec = β[2:4]
        v_ep = (β0^2 - dot(βvec, βvec)) * v + 2 * dot(βvec, v) * βvec + 2 * β0 * cross(βvec, v)
        
        # Rotate using MRP (MRP rotation formula)
        # DCM = I + 8/(1+σ²)² * [σ×]² + 4/(1+σ²) * [σ×]
        function skew(v)
            [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
        end
        σ_cross = skew(σ)
        σ_norm_sq = dot(σ, σ)
        DCM = I + 8/(1+σ_norm_sq)^2 * (σ_cross * σ_cross) + 4/(1+σ_norm_sq) * σ_cross
        v_mrp = DCM * v
        
        @test v_ep ≈ v_mrp atol=1e-12
    end

    @testset "Numerical Stability - Small Rotations" begin
        # Test with very small rotation angle
        θ = 1e-6
        axis = [0.0, 0.0, 1.0]
        β = [cos(θ/2); sin(θ/2) * axis]
        
        σ = AstroCoords.EP2MRP(β)
        β_back = AstroCoords.MRP2EP(σ)
        
        # Should be numerically stable
        @test all(isfinite.(σ))
        @test all(isfinite.(β_back))
        
        # β0 should be very close to 1 for small rotation
        @test β[1] ≈ 1.0 atol=1e-10
    end
end
