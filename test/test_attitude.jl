using Test
using AstroCoords
using LinearAlgebra

@testset "Attitude Transformations" begin
    @testset "EP ↔ MRP conversions" begin
        @testset "Identity rotation" begin
            # β = [1, 0, 0, 0] represents identity (no rotation)
            β = [1.0, 0.0, 0.0, 0.0]
            σ = EP2MRP(β)
            
            # MRP for identity should be zero vector
            @test σ ≈ [0.0, 0.0, 0.0] atol=1e-15
            
            # Round-trip should recover identity
            β_recovered = MRP2EP(σ)
            @test β_recovered ≈ β atol=1e-15
        end
        
        @testset "90° rotation about x-axis" begin
            # β = [cos(45°), sin(45°), 0, 0] for 90° about x
            θ = π/2  # 90 degrees
            β = [cos(θ/2), sin(θ/2), 0.0, 0.0]
            
            σ = EP2MRP(β)
            
            # For 90° rotation: |σ| = tan(θ/4) = tan(π/8) ≈ 0.4142
            σ_expected_mag = tan(θ/4)
            @test norm(σ) ≈ σ_expected_mag atol=1e-10 rtol=1e-10
            
            # Direction should be along x-axis
            @test σ[2] ≈ 0.0 atol=1e-15
            @test σ[3] ≈ 0.0 atol=1e-15
            
            # Round-trip
            β_recovered = MRP2EP(σ)
            @test β_recovered ≈ β atol=1e-12 rtol=1e-10
        end
        
        @testset "90° rotation about y-axis" begin
            θ = π/2
            β = [cos(θ/2), 0.0, sin(θ/2), 0.0]
            
            σ = EP2MRP(β)
            σ_expected_mag = tan(θ/4)
            @test norm(σ) ≈ σ_expected_mag atol=1e-10 rtol=1e-10
            
            @test σ[1] ≈ 0.0 atol=1e-15
            @test σ[3] ≈ 0.0 atol=1e-15
            
            β_recovered = MRP2EP(σ)
            @test β_recovered ≈ β atol=1e-12 rtol=1e-10
        end
        
        @testset "90° rotation about z-axis" begin
            θ = π/2
            β = [cos(θ/2), 0.0, 0.0, sin(θ/2)]
            
            σ = EP2MRP(β)
            σ_expected_mag = tan(θ/4)
            @test norm(σ) ≈ σ_expected_mag atol=1e-10 rtol=1e-10
            
            @test σ[1] ≈ 0.0 atol=1e-15
            @test σ[2] ≈ 0.0 atol=1e-15
            
            β_recovered = MRP2EP(σ)
            @test β_recovered ≈ β atol=1e-12 rtol=1e-10
        end
        
        @testset "180° rotations (β0 = 0)" begin
            # 180° rotation about x-axis
            β_x = [0.0, 1.0, 0.0, 0.0]
            σ_x = EP2MRP(β_x)
            
            # At β0=0, the conversion has special handling
            # |σ| = tan(π/4) = 1 for 180° rotation
            @test norm(σ_x) ≈ 1.0 atol=1e-10
            @test σ_x[2] ≈ 0.0 atol=1e-15
            @test σ_x[3] ≈ 0.0 atol=1e-15
            
            β_x_recovered = MRP2EP(σ_x)
            # Note: EP has sign ambiguity (β and -β represent same rotation)
            @test abs.(β_x_recovered) ≈ abs.(β_x) atol=1e-12
            
            # 180° rotation about y-axis
            β_y = [0.0, 0.0, 1.0, 0.0]
            σ_y = EP2MRP(β_y)
            @test norm(σ_y) ≈ 1.0 atol=1e-10
            @test σ_y[1] ≈ 0.0 atol=1e-15
            @test σ_y[3] ≈ 0.0 atol=1e-15
            
            # 180° rotation about z-axis
            β_z = [0.0, 0.0, 0.0, 1.0]
            σ_z = EP2MRP(β_z)
            @test norm(σ_z) ≈ 1.0 atol=1e-10
            @test σ_z[1] ≈ 0.0 atol=1e-15
            @test σ_z[2] ≈ 0.0 atol=1e-15
        end
        
        @testset "Small rotations (|σ| << 1)" begin
            # Small rotation: 5 degrees about [1,1,1] axis
            θ = 5.0 * π/180
            axis = normalize([1.0, 1.0, 1.0])
            β = [cos(θ/2); sin(θ/2) * axis]
            
            σ = EP2MRP(β)
            
            # For small angles: σ ≈ β_vec / (1 + β0) ≈ β_vec / 2
            @test norm(σ) < 0.1  # Should be small
            
            # Round-trip should be very accurate for small rotations
            β_recovered = MRP2EP(σ)
            @test β_recovered ≈ β atol=1e-14 rtol=1e-12
        end
        
        @testset "Singularity approach: β0 → -1 (shadow set)" begin
            # Test behavior near the singularity where shadow set switching occurs
            # This tests the epsilon handling in EP2MRP
            
            # Rotation slightly less than 360° (β0 slightly > -1)
            θ = 2π - 0.01  # Just before 360°
            β = [cos(θ/2), sin(θ/2), 0.0, 0.0]
            # β0 ≈ -0.99999... (very close to -1)
            
            σ = EP2MRP(β)
            
            # MRP magnitude should be large (approaching infinity as θ→2π)
            # tan((2π - 0.01)/4) ≈ tan(π/2 - 0.0025) ≈ 400
            @test norm(σ) > 100.0  # Should be very large
            
            # Round-trip: should use shadow set formulation
            β_recovered = MRP2EP(σ)
            # Due to shadow set, recovered β might differ by sign
            @test abs.(β_recovered) ≈ abs.(β) atol=1e-8 rtol=1e-6
        end
        
        @testset "Large rotation: 270° about arbitrary axis" begin
            θ = 3π/2  # 270 degrees
            axis = normalize([1.0, 2.0, -1.0])
            β = [cos(θ/2); sin(θ/2) * axis]
            
            σ = EP2MRP(β)
            
            # |σ| = tan(3π/8) ≈ 2.414
            σ_expected_mag = tan(θ/4)
            @test norm(σ) ≈ σ_expected_mag atol=1e-10 rtol=1e-8
            
            # Round-trip
            β_recovered = MRP2EP(σ)
            # For large rotations, may need shadow set
            @test norm(β_recovered - β) < 1e-10 || norm(β_recovered + β) < 1e-10
        end
        
        @testset "Property: EP normalization" begin
            # All Euler parameters should satisfy β0² + β1² + β2² + β3² = 1
            test_rotations = [
                [0.7071, 0.7071, 0.0, 0.0],
                [0.5, 0.5, 0.5, 0.5],
                [0.0, 0.6, 0.8, 0.0],
            ]
            
            for β in test_rotations
                @test sum(β .^ 2) ≈ 1.0 atol=1e-15
                
                σ = EP2MRP(β)
                β_recovered = MRP2EP(σ)
                
                # Recovered EP should also be normalized
                @test sum(β_recovered .^ 2) ≈ 1.0 atol=1e-15
            end
        end
        
        @testset "Property: MRP magnitude relates to rotation angle" begin
            # Test that |σ| = tan(θ/4) for various angles
            test_angles = [π/6, π/4, π/3, π/2, 2π/3, π, 5π/4]
            
            for θ in test_angles
                β = [cos(θ/2), sin(θ/2), 0.0, 0.0]
                σ = EP2MRP(β)
                
                expected_mag = tan(θ/4)
                @test norm(σ) ≈ expected_mag atol=1e-10 rtol=1e-8
            end
        end
    end
end
