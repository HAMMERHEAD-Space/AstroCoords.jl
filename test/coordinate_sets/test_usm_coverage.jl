using Test
using AstroCoords
using StaticArrays

@testset "USM7 Coverage Tests" begin
    @testset "Constructor: AbstractVector{T}" begin
        # Test constructor from generic AbstractVector
        vec = [1.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7071]
        usm = USM7(vec)
        @test usm.C ≈ 1.5 atol=1e-15
        @test usm.Rf1 ≈ 0.1 atol=1e-15
        @test usm.Rf2 ≈ 0.2 atol=1e-15
        @test usm.ϵO1 ≈ 0.3 atol=1e-15
        @test usm.ϵO2 ≈ 0.4 atol=1e-15
        @test usm.ϵO3 ≈ 0.5 atol=1e-15
        @test usm.η0 ≈ 0.7071 atol=1e-15
        @test typeof(usm) === USM7{Float64}
    end

    @testset "Constructor: Individual args" begin
        # Test constructor with individual arguments
        usm = USM7(1.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7071)
        @test usm.C ≈ 1.5 atol=1e-15
        @test usm.Rf1 ≈ 0.1 atol=1e-15
        @test usm.Rf2 ≈ 0.2 atol=1e-15
        @test usm.ϵO1 ≈ 0.3 atol=1e-15
        @test usm.ϵO2 ≈ 0.4 atol=1e-15
        @test usm.ϵO3 ≈ 0.5 atol=1e-15
        @test usm.η0 ≈ 0.7071 atol=1e-15
        
        # Test type promotion
        usm_mixed = USM7(1, 0.1f0, 0.2, 0.3, 0.4, 0.5, 0.7071)
        @test typeof(usm_mixed) === USM7{Float64}
    end

    @testset "Constructor: StaticVector" begin
        # Test constructor from StaticVector
        svec = SVector{7}(1.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7071)
        usm = USM7(svec)
        @test usm.C ≈ 1.5 atol=1e-15
        @test usm.Rf1 ≈ 0.1 atol=1e-15
        @test usm.Rf2 ≈ 0.2 atol=1e-15
        @test usm.ϵO1 ≈ 0.3 atol=1e-15
        @test usm.ϵO2 ≈ 0.4 atol=1e-15
        @test usm.ϵO3 ≈ 0.5 atol=1e-15
        @test usm.η0 ≈ 0.7071 atol=1e-15
        
        # Test USM7{T}(::StaticVector)
        usm_typed = USM7{Float64}(svec)
        @test usm_typed.C ≈ 1.5 atol=1e-15
        @test typeof(usm_typed) === USM7{Float64}
    end

    @testset "params() extraction" begin
        # Test params() returns correct SVector
        usm = USM7(1.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7071)
        p = params(usm)
        @test typeof(p) === SVector{7, Float64}
        @test p[1] ≈ 1.5 atol=1e-15
        @test p[2] ≈ 0.1 atol=1e-15
        @test p[3] ≈ 0.2 atol=1e-15
        @test p[4] ≈ 0.3 atol=1e-15
        @test p[5] ≈ 0.4 atol=1e-15
        @test p[6] ≈ 0.5 atol=1e-15
        @test p[7] ≈ 0.7071 atol=1e-15
    end

    @testset "Edge case: Circular orbit (Rf1=Rf2=0, e=0)" begin
        # Circular orbit: Rf1=Rf2=0, R=0, e=0
        usm = USM7(2.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.9)
        @test usm.C ≈ 2.0 atol=1e-15
        @test usm.Rf1 ≈ 0.0 atol=1e-15
        @test usm.Rf2 ≈ 0.0 atol=1e-15
        
        p = params(usm)
        @test p[2] ≈ 0.0 atol=1e-15
        @test p[3] ≈ 0.0 atol=1e-15
    end

    @testset "Edge case: C→0 limit (low angular momentum)" begin
        # Very small C (approaching zero angular momentum limit)
        usm = USM7(1e-6, 0.01, 0.02, 0.1, 0.2, 0.3, 0.9)
        @test usm.C ≈ 1e-6 atol=1e-15
        @test usm.C > 0.0  # C should be positive
        
        p = params(usm)
        @test p[1] ≈ 1e-6 atol=1e-15
    end

    @testset "Quaternion normalization constraint" begin
        # Test that quaternion components can represent normalized state
        # Note: Constructor doesn't enforce normalization, but we test it accepts values
        ϵO1, ϵO2, ϵO3, η0 = 0.3, 0.4, 0.5, sqrt(1.0 - 0.3^2 - 0.4^2 - 0.5^2)
        usm = USM7(1.5, 0.1, 0.2, ϵO1, ϵO2, ϵO3, η0)
        
        # Verify quaternion magnitude
        quat_mag = sqrt(usm.ϵO1^2 + usm.ϵO2^2 + usm.ϵO3^2 + usm.η0^2)
        @test quat_mag ≈ 1.0 atol=1e-10
    end

    @testset "Type conversion edge cases" begin
        # Test Float32
        usm_f32 = USM7(Float32[1.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7])
        @test typeof(usm_f32) === USM7{Float32}
        
        # Test mixed precision
        usm_mixed = USM7(1.5f0, 0.1, 0.2f0, 0.3, 0.4f0, 0.5, 0.7)
        @test typeof(usm_mixed) === USM7{Float64}
    end
end

@testset "USM6 Coverage Tests" begin
    @testset "Constructor: AbstractVector{T}" begin
        vec = [1.5, 0.1, 0.2, 0.15, 0.25, 0.35]
        usm6 = USM6(vec)
        @test usm6.C ≈ 1.5 atol=1e-15
        @test usm6.Rf1 ≈ 0.1 atol=1e-15
        @test usm6.Rf2 ≈ 0.2 atol=1e-15
        @test usm6.σ1 ≈ 0.15 atol=1e-15
        @test usm6.σ2 ≈ 0.25 atol=1e-15
        @test usm6.σ3 ≈ 0.35 atol=1e-15
        @test typeof(usm6) === USM6{Float64}
    end

    @testset "Constructor: StaticVector" begin
        svec = SVector{6}(1.5, 0.1, 0.2, 0.15, 0.25, 0.35)
        usm6 = USM6(svec)
        @test usm6.C ≈ 1.5 atol=1e-15
        
        usm6_typed = USM6{Float64}(svec)
        @test typeof(usm6_typed) === USM6{Float64}
    end

    @testset "params() extraction" begin
        usm6 = USM6(1.5, 0.1, 0.2, 0.15, 0.25, 0.35)
        p = params(usm6)
        @test typeof(p) === SVector{6, Float64}
        @test p[1] ≈ 1.5 atol=1e-15
        @test p[4] ≈ 0.15 atol=1e-15
    end
end

@testset "USMEM Coverage Tests" begin
    @testset "Constructor: AbstractVector{T}" begin
        vec = [1.5, 0.1, 0.2, 0.05, 0.06, 0.07]
        usmem = USMEM(vec)
        @test usmem.C ≈ 1.5 atol=1e-15
        @test usmem.Rf1 ≈ 0.1 atol=1e-15
        @test usmem.Rf2 ≈ 0.2 atol=1e-15
        @test usmem.a1 ≈ 0.05 atol=1e-15
        @test usmem.a2 ≈ 0.06 atol=1e-15
        @test usmem.a3 ≈ 0.07 atol=1e-15
        @test typeof(usmem) === USMEM{Float64}
    end

    @testset "Constructor: StaticVector" begin
        svec = SVector{6}(1.5, 0.1, 0.2, 0.05, 0.06, 0.07)
        usmem = USMEM(svec)
        @test usmem.C ≈ 1.5 atol=1e-15
        
        usmem_typed = USMEM{Float64}(svec)
        @test typeof(usmem_typed) === USMEM{Float64}
    end

    @testset "params() extraction" begin
        usmem = USMEM(1.5, 0.1, 0.2, 0.05, 0.06, 0.07)
        p = params(usmem)
        @test typeof(p) === SVector{6, Float64}
        @test p[1] ≈ 1.5 atol=1e-15
        @test p[4] ≈ 0.05 atol=1e-15
    end
end
