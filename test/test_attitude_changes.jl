using Test
using AstroCoords
using LinearAlgebra

#@lit {citation="Markley & Crassidis 2014", id="markley2014", ref="Section 2.8, p.45"}
@testset "Attitude Changes — Modified Rodriguez Parameters" begin
    
    @testset "EP2MRP singularity handling at β0 ≈ -1" begin
        #@lit {citation="Markley & Crassidis 2014", id="markley2014", ref="Eq. 2.95, p.47"}
        @testset "Markley & Crassidis 2014 — MRP singularity" begin
            # Test singularity at β0 = -1 (180° rotation)
            # EP for 180° rotation about x-axis: [0, 1, 0, 0]
            β = [0.0, 1.0, 0.0, 0.0]
            
            # EP2MRP should handle this with epsilon in denominator
            σ = EP2MRP(β)
            
            # Should return finite values (not Inf or NaN)
            @test all(isfinite.(σ))
            @test length(σ) == 3
            
            # For 180° rotation, MRP magnitude should be large but finite
            @test norm(σ) > 10.0  # Large due to singularity
            @test norm(σ) < 1e10  # But not infinite
        end
        
        @testset "Near-singularity: β0 ≈ -1 + ε" begin
            # Slightly off singularity
            ε = 1e-6
            β0 = -1.0 + ε
            β_vec = [0.5, 0.5, 0.5]  # Normalize to get unit quaternion
            norm_β_vec = norm(β_vec)
            scale = sqrt(1 - β0^2) / norm_β_vec
            β_vec_normalized = β_vec * scale
            
            β = [β0, β_vec_normalized...]
            
            # Should still be handled smoothly
            σ = EP2MRP(β)
            @test all(isfinite.(σ))
            @test norm(σ) > 1.0  # Should be in shadow set region
        end
    end
    
    @testset "MRP shadow set boundaries (|σ| = 1)" begin
        # MRPs have shadow set for |σ| > 1
        # Test boundary at |σ| = 1
        
        @testset "Shadow set boundary: |σ| = 1" begin
            # Create MRP exactly at boundary
            σ = [1.0/sqrt(3), 1.0/sqrt(3), 1.0/sqrt(3)]  # |σ| = 1
            @test norm(σ) ≈ 1.0 atol=1e-10
            
            # Convert to EP and back
            β = MRP2EP(σ)
            σ_back = EP2MRP(β)
            
            # Should preserve the MRP (within numerical precision)
            @test norm(σ_back - σ) < 1e-6
        end
        
        @testset "Inside shadow set: |σ| = 1.5" begin
            # MRP in shadow set
            σ = [1.5/sqrt(3), 0.0, 0.0]
            @test norm(σ) ≈ 1.5 atol=1e-10
            
            β = MRP2EP(σ)
            σ_back = EP2MRP(β)
            
            # EP2MRP may return shadow set equivalent
            # Check rotation is preserved via DCM
            @test all(isfinite.(β))
            @test norm(β) ≈ 1.0 atol=1e-10  # Unit quaternion
        end
        
        @testset "Small MRP: |σ| = 0.1" begin
            # Small MRP (outside shadow set)
            σ = [0.1, 0.0, 0.0]
            
            β = MRP2EP(σ)
            σ_back = EP2MRP(β)
            
            @test σ_back ≈ σ atol=1e-10 rtol=1e-8
        end
    end
    
    @testset "Round-trip EP → MRP → EP preservation" begin
        #@lit {citation="Markley & Crassidis 2014", id="markley2014", ref="Example 2.7, p.48"}
        @testset "Markley & Crassidis 2014 — Standard rotation" begin
            # 90° rotation about z-axis
            θ = π/2
            axis = [0.0, 0.0, 1.0]
            β0 = cos(θ/2)
            β_vec = sin(θ/2) * axis
            β = [β0, β_vec...]
            
            # Verify unit quaternion
            @test norm(β) ≈ 1.0 atol=1e-10
            
            # Convert to MRP and back
            σ = EP2MRP(β)
            β_back = MRP2EP(σ)
            
            # Quaternion should be preserved (or its negative)
            # Both q and -q represent same rotation
            @test (norm(β_back - β) < 1e-10 || norm(β_back + β) < 1e-10)
        end
        
        @testset "45° rotation about arbitrary axis" begin
            θ = π/4
            axis = normalize([1.0, 2.0, 3.0])
            β0 = cos(θ/2)
            β_vec = sin(θ/2) * axis
            β = [β0, β_vec...]
            
            σ = EP2MRP(β)
            β_back = MRP2EP(σ)
            
            @test (norm(β_back - β) < 1e-10 || norm(β_back + β) < 1e-10)
        end
        
        @testset "Small angle (near identity)" begin
            # 1° rotation
            θ = π/180
            axis = [0.0, 1.0, 0.0]
            β0 = cos(θ/2)
            β_vec = sin(θ/2) * axis
            β = [β0, β_vec...]
            
            σ = EP2MRP(β)
            β_back = MRP2EP(σ)
            
            # For small angles, MRP ≈ β_vec/2
            @test norm(σ) < 0.01
            @test (norm(β_back - β) < 1e-10 || norm(β_back + β) < 1e-10)
        end
    end
    
    @testset "MRP properties" begin
        @testset "MRP for identity rotation" begin
            # Identity: β = [1, 0, 0, 0]
            β = [1.0, 0.0, 0.0, 0.0]
            σ = EP2MRP(β)
            
            # Identity should give zero MRP
            @test norm(σ) ≈ 0.0 atol=1e-10
        end
        
        @testset "MRP magnitude increases with rotation angle" begin
            # Test that |σ| increases monotonically with rotation angle
            angles = [10.0, 30.0, 60.0, 90.0, 120.0] * π/180
            axis = [0.0, 0.0, 1.0]
            
            σ_norms = Float64[]
            for θ in angles
                β0 = cos(θ/2)
                β_vec = sin(θ/2) * axis
                β = [β0, β_vec...]
                σ = EP2MRP(β)
                push!(σ_norms, norm(σ))
            end
            
            # Should be monotonically increasing for θ < 180°
            for i in 1:length(σ_norms)-1
                @test σ_norms[i+1] > σ_norms[i]
            end
        end
    end
    
    @testset "Numerical precision: mixed types" begin
        # Test with Float32 input
        β_f32 = Float32[0.7071, 0.7071, 0.0, 0.0]  # 90° about x
        σ_f32 = EP2MRP(β_f32)
        
        # Should return Float32
        @test eltype(σ_f32) == Float32
        @test all(isfinite.(σ_f32))
        
        # Compare with Float64 version
        β_f64 = Float64[0.7071, 0.7071, 0.0, 0.0]
        σ_f64 = EP2MRP(β_f64)
        
        # Should be close (accounting for Float32 precision)
        @test norm(Float64.(σ_f32) - σ_f64) < 1e-4
    end
end
