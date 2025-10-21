using Test
using AstroCoords
using StaticArrays

@testset "KustaanheimoStiefel Coverage Tests" begin
    @testset "Constructor: AbstractVector{T}" begin
        # Test constructor from generic AbstractVector (currently 0% covered)
        vec = [1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, -0.5, 100.0]
        ks = KustaanheimoStiefel(vec)
        @test ks.u₁ ≈ 1.0 atol=1e-15
        @test ks.u₂ ≈ 2.0 atol=1e-15
        @test ks.u₃ ≈ 3.0 atol=1e-15
        @test ks.u₄ ≈ 4.0 atol=1e-15
        @test ks.u₁_prime ≈ 0.5 atol=1e-15
        @test ks.u₂_prime ≈ 0.6 atol=1e-15
        @test ks.u₃_prime ≈ 0.7 atol=1e-15
        @test ks.u₄_prime ≈ 0.8 atol=1e-15
        @test ks.h ≈ -0.5 atol=1e-15
        @test ks.τ ≈ 100.0 atol=1e-15
        @test typeof(ks) === KustaanheimoStiefel{Float64}
    end

    @testset "Constructor: Individual args" begin
        # Test constructor with individual arguments (currently low coverage)
        ks = KustaanheimoStiefel(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, -0.5, 100.0)
        @test ks.u₁ ≈ 1.0 atol=1e-15
        @test ks.u₂ ≈ 2.0 atol=1e-15
        @test ks.u₃ ≈ 3.0 atol=1e-15
        @test ks.u₄ ≈ 4.0 atol=1e-15
        @test ks.u₁_prime ≈ 0.5 atol=1e-15
        @test ks.u₂_prime ≈ 0.6 atol=1e-15
        @test ks.u₃_prime ≈ 0.7 atol=1e-15
        @test ks.u₄_prime ≈ 0.8 atol=1e-15
        @test ks.h ≈ -0.5 atol=1e-15
        @test ks.τ ≈ 100.0 atol=1e-15
        @test typeof(ks) === KustaanheimoStiefel{Float64}
        
        # Test type promotion with mixed types
        ks_mixed = KustaanheimoStiefel(1, 2.0f0, 3.0, 4, 0.5, 0.6f0, 0.7, 0.8, -0.5, 100.0)
        @test typeof(ks_mixed) === KustaanheimoStiefel{Float64}
    end

    @testset "Constructor: StaticVector" begin
        # Test KustaanheimoStiefel(::StaticVector) constructor
        svec = SVector{10}(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, -0.5, 100.0)
        ks = KustaanheimoStiefel(svec)
        @test ks.u₁ ≈ 1.0 atol=1e-15
        @test ks.u₂ ≈ 2.0 atol=1e-15
        @test typeof(ks) === KustaanheimoStiefel{Float64}
        
        # Test KustaanheimoStiefel{T}(::StaticVector) with explicit type
        ks_typed = KustaanheimoStiefel{Float64}(svec)
        @test ks_typed.u₁ ≈ 1.0 atol=1e-15
        @test typeof(ks_typed) === KustaanheimoStiefel{Float64}
        
        # Test with Float32 StaticVector
        svec_f32 = SVector{10,Float32}(1.0f0, 2.0f0, 3.0f0, 4.0f0, 0.5f0, 0.6f0, 0.7f0, 0.8f0, -0.5f0, 100.0f0)
        ks_f32 = KustaanheimoStiefel(svec_f32)
        @test typeof(ks_f32) === KustaanheimoStiefel{Float32}
    end

    @testset "params() extraction" begin
        # Test params() returns correct SVector
        ks = KustaanheimoStiefel(1.1, 2.2, 3.3, 4.4, 0.55, 0.66, 0.77, 0.88, -0.55, 110.0)
        p = params(ks)
        @test typeof(p) === SVector{10, Float64}
        @test p[1] ≈ 1.1 atol=1e-15
        @test p[2] ≈ 2.2 atol=1e-15
        @test p[3] ≈ 3.3 atol=1e-15
        @test p[4] ≈ 4.4 atol=1e-15
        @test p[5] ≈ 0.55 atol=1e-15
        @test p[6] ≈ 0.66 atol=1e-15
        @test p[7] ≈ 0.77 atol=1e-15
        @test p[8] ≈ 0.88 atol=1e-15
        @test p[9] ≈ -0.55 atol=1e-15
        @test p[10] ≈ 110.0 atol=1e-15
    end

    @testset "Edge case: Zero energy (h=0)" begin
        # Test with h = 0 (parabolic orbit at energy boundary)
        ks = KustaanheimoStiefel(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, 0.0, 100.0)
        @test ks.h ≈ 0.0 atol=1e-15
        
        # Verify params extraction works
        p = params(ks)
        @test p[9] ≈ 0.0 atol=1e-15
    end

    @testset "Edge case: Negative energy (h<0, bound orbit)" begin
        # Test with h < 0 (bound/elliptic orbit)
        ks = KustaanheimoStiefel(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, -1.5, 100.0)
        @test ks.h ≈ -1.5 atol=1e-15
        @test ks.h < 0.0
    end

    @testset "Edge case: Positive energy (h>0, unbound orbit)" begin
        # Test with h > 0 (unbound/hyperbolic orbit)
        ks = KustaanheimoStiefel(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, 0.8, 100.0)
        @test ks.h ≈ 0.8 atol=1e-15
        @test ks.h > 0.0
    end

    @testset "Edge case: τ initialization (time element)" begin
        # Test with τ = 0 (initial time at origin)
        ks = KustaanheimoStiefel(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, -0.5, 0.0)
        @test ks.τ ≈ 0.0 atol=1e-15
        
        # Test with negative τ
        ks_neg = KustaanheimoStiefel(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, -0.5, -50.0)
        @test ks_neg.τ ≈ -50.0 atol=1e-15
    end

    @testset "Edge case: u-vector normalization constraints" begin
        # Test KS position vector (u₁, u₂, u₃, u₄) with specific magnitudes
        # For KS transformation: r = u₁² + u₂² + u₃² + u₄²
        u_mag_sq = 1.0^2 + 2.0^2 + 3.0^2 + 4.0^2  # = 30.0
        ks = KustaanheimoStiefel(1.0, 2.0, 3.0, 4.0, 0.5, 0.6, 0.7, 0.8, -0.5, 100.0)
        
        # Calculate position magnitude from KS components
        r_from_ks = ks.u₁^2 + ks.u₂^2 + ks.u₃^2 + ks.u₄^2
        @test r_from_ks ≈ 30.0 atol=1e-12
    end

    @testset "Edge case: Small u-vector components" begin
        # Test with very small u-vector components (near-singularity)
        ks = KustaanheimoStiefel(1e-6, 1e-6, 1e-6, 1e-6, 0.1, 0.1, 0.1, 0.1, -0.5, 100.0)
        @test ks.u₁ ≈ 1e-6 atol=1e-15
        @test ks.u₂ ≈ 1e-6 atol=1e-15
        
        # Position magnitude should be very small
        r_from_ks = ks.u₁^2 + ks.u₂^2 + ks.u₃^2 + ks.u₄^2
        @test r_from_ks ≈ 4e-12 atol=1e-15
    end

    @testset "Type stability checks" begin
        # Verify type stability across different input types
        vec_int = [1, 2, 3, 4, 0, 0, 0, 0, 0, 0]
        ks_int = KustaanheimoStiefel(vec_int)
        @test typeof(ks_int) === KustaanheimoStiefel{Int64}
        
        vec_f64 = [1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ks_f64 = KustaanheimoStiefel(vec_f64)
        @test typeof(ks_f64) === KustaanheimoStiefel{Float64}
    end
end
