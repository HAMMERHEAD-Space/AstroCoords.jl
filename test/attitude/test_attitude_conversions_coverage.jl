using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "Attitude Conversions Coverage" begin
    @testset "EP2MRP - identity rotation" begin
        # Identity rotation: β = [1, 0, 0, 0]
        β_identity = [1.0, 0.0, 0.0, 0.0]
        σ = EP2MRP(β_identity)
        
        # MRP for identity should be zero vector
        @test σ[1] ≈ 0.0 atol=1e-14
        @test σ[2] ≈ 0.0 atol=1e-14
        @test σ[3] ≈ 0.0 atol=1e-14
        
        # Check return type
        @test σ isa SVector{3}
    end
    
    @testset "EP2MRP - 90° rotation about X-axis" begin
        # 90° about X: β = [cos(45°), sin(45°), 0, 0]
        β_x = [cos(π/4), sin(π/4), 0.0, 0.0]
        σ = EP2MRP(β_x)
        
        # For 90° rotation, MRP components follow tan(θ/4)*axis
        # θ=π/2, so tan(π/8) ≈ 0.4142
        expected_mag = tan(π/8)
        @test σ[1] ≈ expected_mag atol=1e-10
        @test σ[2] ≈ 0.0 atol=1e-14
        @test σ[3] ≈ 0.0 atol=1e-14
    end
    
    @testset "EP2MRP - 90° rotation about Y-axis" begin
        β_y = [cos(π/4), 0.0, sin(π/4), 0.0]
        σ = EP2MRP(β_y)
        
        expected_mag = tan(π/8)
        @test σ[1] ≈ 0.0 atol=1e-14
        @test σ[2] ≈ expected_mag atol=1e-10
        @test σ[3] ≈ 0.0 atol=1e-14
    end
    
    @testset "EP2MRP - 90° rotation about Z-axis" begin
        β_z = [cos(π/4), 0.0, 0.0, sin(π/4)]
        σ = EP2MRP(β_z)
        
        expected_mag = tan(π/8)
        @test σ[1] ≈ 0.0 atol=1e-14
        @test σ[2] ≈ 0.0 atol=1e-14
        @test σ[3] ≈ expected_mag atol=1e-10
    end
    
    @testset "EP2MRP - 180° rotation (near singularity)" begin
        # 180° rotation about X: β = [0, 1, 0, 0]
        β_180 = [0.0, 1.0, 0.0, 0.0]
        σ = EP2MRP(β_180)
        
        # At 180°, MRPs become infinite, but we should get large values
        @test abs(σ[1]) > 1.0  # Should be large
        @test !isinf(σ[1])      # Should not be infinite (numerical stability)
        @test !isnan(σ[1])      # Should not be NaN
    end
    
    @testset "EP2MRP - near-singularity (β0 ≈ -1)" begin
        # β0 close to -1 (near 360° rotation / shadow set)
        β_near = [-0.99, 0.1, 0.1, 0.1]
        β_near = β_near ./ norm(β_near)  # Normalize
        σ = EP2MRP(β_near)
        
        # Should handle gracefully without NaN/Inf
        @test !any(isnan.(σ))
        @test !any(isinf.(σ))
        @test norm(σ) > 1.0  # Should be large near singularity
    end
    
    @testset "EP2MRP - small rotations" begin
        # Small rotation (5°) about arbitrary axis
        θ = 5 * π/180
        axis = normalize([1.0, 1.0, 1.0])
        β_small = [cos(θ/2), axis[1]*sin(θ/2), axis[2]*sin(θ/2), axis[3]*sin(θ/2)]
        σ = EP2MRP(β_small)
        
        # For small angles, MRP ≈ (1/2)*β_vector
        expected_approx = 0.5 .* β_small[2:4]
        @test σ[1] ≈ expected_approx[1] atol=1e-3
        @test σ[2] ≈ expected_approx[2] atol=1e-3
        @test σ[3] ≈ expected_approx[3] atol=1e-3
    end
    
    @testset "MRP2EP - zero MRP (identity rotation)" begin
        σ_zero = [0.0, 0.0, 0.0]
        β = MRP2EP(σ_zero)
        
        # Should give identity EP
        @test β[1] ≈ 1.0 atol=1e-14
        @test β[2] ≈ 0.0 atol=1e-14
        @test β[3] ≈ 0.0 atol=1e-14
        @test β[4] ≈ 0.0 atol=1e-14
        
        # Check unit quaternion constraint
        @test norm(β) ≈ 1.0 atol=1e-14
    end
    
    @testset "MRP2EP - 90° rotations" begin
        # Approximate MRP for 90° about X
        σ_x = [tan(π/8), 0.0, 0.0]
        β = MRP2EP(σ_x)
        
        # Should be close to [cos(π/4), sin(π/4), 0, 0]
        @test β[1] ≈ cos(π/4) atol=1e-10
        @test β[2] ≈ sin(π/4) atol=1e-10
        @test β[3] ≈ 0.0 atol=1e-14
        @test β[4] ≈ 0.0 atol=1e-14
        @test norm(β) ≈ 1.0 atol=1e-12
    end
    
    @testset "MRP2EP - shadow set handling" begin
        # Large MRP (beyond singularity, should use shadow set)
        σ_large = [2.0, 0.0, 0.0]
        β = MRP2EP(σ_large)
        
        # Should still be valid unit quaternion
        @test norm(β) ≈ 1.0 atol=1e-12
        @test !any(isnan.(β))
        @test !any(isinf.(β))
    end
    
    @testset "Round-trip: EP → MRP → EP" begin
        # Test various rotations
        test_cases = [
            [1.0, 0.0, 0.0, 0.0],  # Identity
            [cos(π/6), sin(π/6), 0.0, 0.0],  # 60° about X
            [cos(π/3), 0.0, sin(π/3), 0.0],  # 120° about Y
            [cos(π/4), 0.0, 0.0, sin(π/4)],  # 90° about Z
            [0.5, 0.5, 0.5, 0.5],  # General rotation
        ]
        
        for β_orig in test_cases
            β_norm = β_orig ./ norm(β_orig)
            σ = EP2MRP(β_norm)
            β_back = MRP2EP(σ)
            
            # Should recover original (or negative, since ±β represent same rotation)
            matches_positive = all(abs.(β_back .- β_norm) .< 1e-10)
            matches_negative = all(abs.(β_back .+ β_norm) .< 1e-10)
            @test matches_positive || matches_negative
        end
    end
    
    @testset "Round-trip: MRP → EP → MRP" begin
        # Test various MRPs
        test_mrps = [
            [0.0, 0.0, 0.0],      # Identity
            [0.1, 0.0, 0.0],      # Small X
            [0.0, 0.2, 0.0],      # Small Y
            [0.0, 0.0, 0.3],      # Small Z
            [0.1, 0.2, 0.3],      # General small
            [0.5, 0.5, 0.5],      # Moderate
        ]
        
        for σ_orig in test_mrps
            β = MRP2EP(σ_orig)
            σ_back = EP2MRP(β)
            
            # Should recover original MRP
            @test σ_back[1] ≈ σ_orig[1] atol=1e-10
            @test σ_back[2] ≈ σ_orig[2] atol=1e-10
            @test σ_back[3] ≈ σ_orig[3] atol=1e-10
        end
    end
    
    @testset "Type stability and promotion" begin
        # Float32 input
        β_f32 = Float32[1.0, 0.0, 0.0, 0.0]
        σ_f32 = EP2MRP(β_f32)
        @test eltype(σ_f32) == Float32
        
        # Mixed types
        β_mixed = [Float64(1.0), Float32(0.0), Float64(0.0), Float32(0.0)]
        σ_mixed = EP2MRP(β_mixed)
        @test σ_mixed isa SVector{3}
        
        # SVector input
        β_sv = SVector{4}(1.0, 0.0, 0.0, 0.0)
        σ_sv = EP2MRP(β_sv)
        @test σ_sv isa SVector{3}
    end
    
    @testset "MRP unit quaternion constraint" begin
        # All conversions should preserve unit norm
        σ_test = [0.3, 0.4, 0.5]
        β = MRP2EP(σ_test)
        @test norm(β) ≈ 1.0 atol=1e-12
        
        # Try with larger MRPs
        σ_large = [1.0, 1.0, 1.0]
        β_large = MRP2EP(σ_large)
        @test norm(β_large) ≈ 1.0 atol=1e-12
    end
end
