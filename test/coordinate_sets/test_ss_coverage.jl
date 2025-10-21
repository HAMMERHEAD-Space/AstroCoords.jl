using Test
using AstroCoords
using StaticArrays

@testset "StiefelScheifele Coverage Tests" begin
    @testset "Constructor: AbstractVector{T}" begin
        # Test constructor from generic AbstractVector (currently 0% covered)
        vec = [1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, 1.5, 100.0]
        ss = StiefelScheifele(vec)
        @test ss.α1 ≈ 1.0 atol=1e-15
        @test ss.α2 ≈ 2.0 atol=1e-15
        @test ss.α3 ≈ 3.0 atol=1e-15
        @test ss.α4 ≈ 4.0 atol=1e-15
        @test ss.β1 ≈ 0.5 atol=1e-15
        @test ss.β2 ≈ 0.6 atol=1e-15
        @test ss.β3 ≈ 0.7 atol=1e-15
        @test ss.β4 ≈ 0.8 atol=1e-15
        @test ss.ω ≈ 1.5 atol=1e-15
        @test ss.t ≈ 100.0 atol=1e-15
        @test typeof(ss) === StiefelScheifele{Float64}
    end

    @testset "Constructor: Individual args" begin
        # Test constructor with individual arguments (currently low coverage)
        ss = StiefelScheifele(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, 1.5, 100.0)
        @test ss.α1 ≈ 1.0 atol=1e-15
        @test ss.α2 ≈ 2.0 atol=1e-15
        @test ss.α3 ≈ 3.0 atol=1e-15
        @test ss.α4 ≈ 4.0 atol=1e-15
        @test ss.β1 ≈ 0.5 atol=1e-15
        @test ss.β2 ≈ 0.6 atol=1e-15
        @test ss.β3 ≈ 0.7 atol=1e-15
        @test ss.β4 ≈ 0.8 atol=1e-15
        @test ss.ω ≈ 1.5 atol=1e-15
        @test ss.t ≈ 100.0 atol=1e-15
        @test typeof(ss) === StiefelScheifele{Float64}
        
        # Test type promotion with mixed types
        ss_mixed = StiefelScheifele(1, 2.0f0, 3.0, 4, 0.5, 0.6f0, 0.7, 0.8, 1.5, 100.0)
        @test typeof(ss_mixed) === StiefelScheifele{Float64}
    end

    @testset "Constructor: StaticVector" begin
        # Test StiefelScheifele(::StaticVector) constructor
        svec = SVector{10}(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, 1.5, 100.0)
        ss = StiefelScheifele(svec)
        @test ss.α1 ≈ 1.0 atol=1e-15
        @test ss.α2 ≈ 2.0 atol=1e-15
        @test typeof(ss) === StiefelScheifele{Float64}
        
        # Test StiefelScheifele{T}(::StaticVector) with explicit type
        ss_typed = StiefelScheifele{Float64}(svec)
        @test ss_typed.α1 ≈ 1.0 atol=1e-15
        @test typeof(ss_typed) === StiefelScheifele{Float64}
        
        # Test with Float32 StaticVector
        svec_f32 = SVector{10,Float32}(1.0f0, 2.0f0, 3.0f0, 4.0f0, 0.5f0, 0.6f0, 0.7f0, 0.8f0, 1.5f0, 100.0f0)
        ss_f32 = StiefelScheifele(svec_f32)
        @test typeof(ss_f32) === StiefelScheifele{Float32}
    end

    @testset "params() extraction" begin
        # Test params() returns correct SVector
        ss = StiefelScheifele(1.1, 2.2, 3.3, 4.4, 0.55, 0.66, 0.77, 0.88, 1.65, 110.0)
        p = params(ss)
        @test typeof(p) === SVector{10, Float64}
        @test p[1] ≈ 1.1 atol=1e-15
        @test p[2] ≈ 2.2 atol=1e-15
        @test p[3] ≈ 3.3 atol=1e-15
        @test p[4] ≈ 4.4 atol=1e-15
        @test p[5] ≈ 0.55 atol=1e-15
        @test p[6] ≈ 0.66 atol=1e-15
        @test p[7] ≈ 0.77 atol=1e-15
        @test p[8] ≈ 0.88 atol=1e-15
        @test p[9] ≈ 1.65 atol=1e-15
        @test p[10] ≈ 110.0 atol=1e-15
    end

    @testset "Edge case: α quaternion pair" begin
        # Test α components (α1, α2, α3, α4) representing quaternion-like state
        # Verify they can represent normalized state
        α_norm = sqrt(1.0^2 + 2.0^2 + 3.0^2 + 4.0^2)
        ss = StiefelScheifele(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, 1.5, 100.0)
        
        α_mag = sqrt(ss.α1^2 + ss.α2^2 + ss.α3^2 + ss.α4^2)
        @test α_mag ≈ α_norm atol=1e-12
    end

    @testset "Edge case: β quaternion pair" begin
        # Test β components (β1, β2, β3, β4) representing quaternion-like state
        β_norm = sqrt(0.5^2 + 0.6^2 + 0.7^2 + 0.8^2)
        ss = StiefelScheifele(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, 1.5, 100.0)
        
        β_mag = sqrt(ss.β1^2 + ss.β2^2 + ss.β3^2 + ss.β4^2)
        @test β_mag ≈ β_norm atol=1e-12
    end

    @testset "Edge case: ω parameter (energy-related)" begin
        # Test ω parameter with various values
        # ω > 0 (positive energy)
        ss_pos = StiefelScheifele(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, 2.0, 100.0)
        @test ss_pos.ω ≈ 2.0 atol=1e-15
        @test ss_pos.ω > 0.0
        
        # ω = 0 (zero energy boundary)
        ss_zero = StiefelScheifele(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, 0.0, 100.0)
        @test ss_zero.ω ≈ 0.0 atol=1e-15
        
        # Very small ω (near-zero energy)
        ss_small = StiefelScheifele(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, 1e-10, 100.0)
        @test ss_small.ω ≈ 1e-10 atol=1e-15
    end

    @testset "Edge case: Time coordinate t handling" begin
        # Test with t = 0 (initial time at origin)
        ss_t0 = StiefelScheifele(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, 1.5, 0.0)
        @test ss_t0.t ≈ 0.0 atol=1e-15
        
        # Test with negative t
        ss_neg = StiefelScheifele(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, 1.5, -50.0)
        @test ss_neg.t ≈ -50.0 atol=1e-15
        
        # Test with large t
        ss_large = StiefelScheifele(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, 1.5, 1e6)
        @test ss_large.t ≈ 1e6 atol=1e-9
    end

    @testset "Edge case: Zero α components" begin
        # Test with zero α components (singular case)
        ss = StiefelScheifele(0.0, 0.0, 0.0, 0.0, 0.5, 0.6, 0.7, 0.8, 1.5, 100.0)
        @test ss.α1 ≈ 0.0 atol=1e-15
        @test ss.α2 ≈ 0.0 atol=1e-15
        @test ss.α3 ≈ 0.0 atol=1e-15
        @test ss.α4 ≈ 0.0 atol=1e-15
    end

    @testset "Edge case: Zero β components" begin
        # Test with zero β components (singular case)
        ss = StiefelScheifele(1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.5, 100.0)
        @test ss.β1 ≈ 0.0 atol=1e-15
        @test ss.β2 ≈ 0.0 atol=1e-15
        @test ss.β3 ≈ 0.0 atol=1e-15
        @test ss.β4 ≈ 0.0 atol=1e-15
    end

    @testset "Type stability checks" begin
        # Verify type stability across different input types
        vec_int = [1, 2, 3, 4, 0, 0, 0, 0, 1, 0]
        ss_int = StiefelScheifele(vec_int)
        @test typeof(ss_int) === StiefelScheifele{Int64}
        
        vec_f64 = [1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        ss_f64 = StiefelScheifele(vec_f64)
        @test typeof(ss_f64) === StiefelScheifele{Float64}
    end

    @testset "Relationship between α and β" begin
        # Test that α and β can represent related quaternion states
        # (They evolve together during integration)
        ss = StiefelScheifele(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0)
        
        # Check orthogonality-like properties
        α_vec = [ss.α1, ss.α2, ss.α3, ss.α4]
        β_vec = [ss.β1, ss.β2, ss.β3, ss.β4]
        
        # Both should have well-defined magnitudes
        α_mag = sqrt(sum(α_vec.^2))
        β_mag = sqrt(sum(β_vec.^2))
        @test α_mag ≈ 1.0 atol=1e-12
        @test β_mag ≈ 1.0 atol=1e-12
    end
end
