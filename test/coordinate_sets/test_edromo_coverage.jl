using Test
using AstroCoords
using StaticArrays

@testset "EDromo Coverage Tests" begin
    @testset "Constructor: AbstractVector{T}" begin
        # Test constructor from generic AbstractVector (currently 0% covered)
        vec = [1.0, 0.5, -0.2, 0.1, 0.2, 0.3, 0.9, 100.0]
        edromo = EDromo(vec)
        @test edromo.ζ₁ ≈ 1.0 atol=1e-15
        @test edromo.ζ₂ ≈ 0.5 atol=1e-15
        @test edromo.ζ₃ ≈ -0.2 atol=1e-15
        @test edromo.ζ₄ ≈ 0.1 atol=1e-15
        @test edromo.ζ₅ ≈ 0.2 atol=1e-15
        @test edromo.ζ₆ ≈ 0.3 atol=1e-15
        @test edromo.ζ₇ ≈ 0.9 atol=1e-15
        @test edromo.ζ₈ ≈ 100.0 atol=1e-15
        @test typeof(edromo) === EDromo{Float64}
    end

    @testset "Constructor: Individual args with type promotion" begin
        # Test constructor with individual arguments (currently low coverage)
        edromo = EDromo(1.0, 0.5, -0.2, 0.1, 0.2, 0.3, 0.9, 100.0)
        @test edromo.ζ₁ ≈ 1.0 atol=1e-15
        @test edromo.ζ₂ ≈ 0.5 atol=1e-15
        @test edromo.ζ₃ ≈ -0.2 atol=1e-15
        @test edromo.ζ₄ ≈ 0.1 atol=1e-15
        @test edromo.ζ₅ ≈ 0.2 atol=1e-15
        @test edromo.ζ₆ ≈ 0.3 atol=1e-15
        @test edromo.ζ₇ ≈ 0.9 atol=1e-15
        @test edromo.ζ₈ ≈ 100.0 atol=1e-15
        @test typeof(edromo) === EDromo{Float64}
        
        # Test type promotion with mixed types (Int, Float32, Float64)
        edromo_mixed = EDromo(1, 0.5f0, -0.2, 0.1, 0.2f0, 0.3, 0.9, 100.0)
        @test typeof(edromo_mixed) === EDromo{Float64}
        @test edromo_mixed.ζ₁ ≈ 1.0 atol=1e-15
    end

    @testset "Constructor: StaticVector edge cases" begin
        # Test EDromo(::StaticVector) constructor
        svec = SVector{8}(1.0, 0.5, -0.2, 0.1, 0.2, 0.3, 0.9, 100.0)
        edromo = EDromo(svec)
        @test edromo.ζ₁ ≈ 1.0 atol=1e-15
        @test edromo.ζ₂ ≈ 0.5 atol=1e-15
        @test typeof(edromo) === EDromo{Float64}
        
        # Test EDromo{T}(::StaticVector) constructor with explicit type
        edromo_typed = EDromo{Float64}(svec)
        @test edromo_typed.ζ₁ ≈ 1.0 atol=1e-15
        @test typeof(edromo_typed) === EDromo{Float64}
        
        # Test with Float32 StaticVector
        svec_f32 = SVector{8,Float32}(1.0f0, 0.5f0, -0.2f0, 0.1f0, 0.2f0, 0.3f0, 0.9f0, 100.0f0)
        edromo_f32 = EDromo(svec_f32)
        @test typeof(edromo_f32) === EDromo{Float32}
    end

    @testset "params() extraction" begin
        # Test params() returns correct SVector
        edromo = EDromo(1.1, 0.6, -0.3, 0.15, 0.25, 0.35, 0.85, 150.0)
        p = params(edromo)
        @test typeof(p) === SVector{8, Float64}
        @test p[1] ≈ 1.1 atol=1e-15
        @test p[2] ≈ 0.6 atol=1e-15
        @test p[3] ≈ -0.3 atol=1e-15
        @test p[4] ≈ 0.15 atol=1e-15
        @test p[5] ≈ 0.25 atol=1e-15
        @test p[6] ≈ 0.35 atol=1e-15
        @test p[7] ≈ 0.85 atol=1e-15
        @test p[8] ≈ 150.0 atol=1e-15
    end

    @testset "Edge case: Quaternion normalization" begin
        # Test quaternion components that should be normalized (ζ₄, ζ₅, ζ₆, ζ₇)
        # Create normalized quaternion components
        q_sum_sq = 0.1^2 + 0.2^2 + 0.3^2
        ζ₇_normalized = sqrt(1.0 - q_sum_sq)
        edromo = EDromo(1.0, 0.5, -0.2, 0.1, 0.2, 0.3, ζ₇_normalized, 100.0)
        
        # Verify quaternion magnitude
        quat_mag = sqrt(edromo.ζ₄^2 + edromo.ζ₅^2 + edromo.ζ₆^2 + edromo.ζ₇^2)
        @test quat_mag ≈ 1.0 atol=1e-10
    end

    @testset "Edge case: Negative energy (ζ₃ < 0, elliptic orbit)" begin
        # Test with negative ζ₃ (corresponds to bound/elliptic orbit)
        edromo = EDromo(1.0, 0.5, -1.5, 0.1, 0.2, 0.3, 0.9, 100.0)
        @test edromo.ζ₃ ≈ -1.5 atol=1e-15
        @test edromo.ζ₃ < 0.0
    end

    @testset "Edge case: Positive energy (ζ₃ > 0, hyperbolic orbit)" begin
        # Test with positive ζ₃ (corresponds to unbound/hyperbolic orbit)
        edromo = EDromo(1.0, 0.5, 0.8, 0.1, 0.2, 0.3, 0.9, 100.0)
        @test edromo.ζ₃ ≈ 0.8 atol=1e-15
        @test edromo.ζ₃ > 0.0
    end

    @testset "Edge case: Zero time element" begin
        # Test with ζ₈ = 0 (initial time at origin)
        edromo = EDromo(1.0, 0.5, -0.2, 0.1, 0.2, 0.3, 0.9, 0.0)
        @test edromo.ζ₈ ≈ 0.0 atol=1e-15
    end

    @testset "Type stability checks" begin
        # Verify type stability across different input types
        vec_int = [1, 0, 0, 0, 0, 0, 1, 0]
        edromo_int = EDromo(vec_int)
        @test typeof(edromo_int) === EDromo{Int64}
        
        vec_f64 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        edromo_f64 = EDromo(vec_f64)
        @test typeof(edromo_f64) === EDromo{Float64}
    end
end
