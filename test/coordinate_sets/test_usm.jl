using Test
using AstroCoords
using StaticArrays

@testset "USM Coordinate Sets" begin
    @testset "USM7 - Inner Constructor" begin
        # Test USM7 inner constructor
        usm7 = USM7{Float64}(1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
        @test usm7 isa USM7{Float64}
        @test usm7.C == 1.0
        @test usm7.Rf1 == 0.1
        @test usm7.Rf2 == 0.2
        @test usm7.ϵO1 == 0.3
        @test usm7.ϵO2 == 0.4
        @test usm7.ϵO3 == 0.5
        @test usm7.η0 == 0.6
    end
    
    @testset "USM7 - Constructor from Vector" begin
        vec = [1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        usm7 = USM7(vec)
        @test usm7 isa USM7{Float64}
        @test usm7.C == 1.0
        @test usm7.Rf1 == 0.1
        @test usm7.Rf2 == 0.2
        @test usm7.ϵO1 == 0.3
        @test usm7.ϵO2 == 0.4
        @test usm7.ϵO3 == 0.5
        @test usm7.η0 == 0.6
        
        # Test with Int vector
        vec_int = [1, 0, 0, 0, 0, 0, 1]
        usm7_int = USM7(vec_int)
        @test usm7_int isa USM7{Int}
    end
    
    @testset "USM7 - Type Promotion" begin
        # Mixed types should promote
        usm7 = USM7(1.0, 0.1f0, 0.2, 0, 0.4, 0.5, 0.6)
        @test usm7 isa USM7
        @test eltype(usm7) <: AbstractFloat
        
        # All Float32
        usm7_f32 = USM7(1.0f0, 0.1f0, 0.2f0, 0.3f0, 0.4f0, 0.5f0, 0.6f0)
        @test usm7_f32 isa USM7{Float32}
    end
    
    @testset "USM7 - StaticVector Constructor" begin
        svec = SVector{7}(1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
        usm7 = USM7(svec)
        @test usm7 isa USM7{Float64}
        @test usm7.C == 1.0
        @test usm7.η0 == 0.6
    end
    
    @testset "USM7 - params() Method" begin
        usm7 = USM7(1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
        p = params(usm7)
        
        @test p isa SVector{7, Float64}
        @test p[1] == 1.0
        @test p[2] == 0.1
        @test p[3] == 0.2
        @test p[4] == 0.3
        @test p[5] == 0.4
        @test p[6] == 0.5
        @test p[7] == 0.6
    end
    
    @testset "USM7 - Base.one()" begin
        usm7_one = one(USM7)
        @test usm7_one isa USM7{Float64}
        @test all(params(usm7_one) .== 0.0)
        
        usm7_one_f32 = one(USM7, T=Float32)
        @test usm7_one_f32 isa USM7{Float32}
    end
    
    @testset "USM7 - Base.getindex" begin
        usm7 = USM7(1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
        
        @test usm7[1] == 1.0
        @test usm7[2] == 0.1
        @test usm7[3] == 0.2
        @test usm7[4] == 0.3
        @test usm7[5] == 0.4
        @test usm7[6] == 0.5
        @test usm7[7] == 0.6
        
        @test_throws BoundsError usm7[0]
        @test_throws BoundsError usm7[8]
    end
    
    @testset "USM6 - Inner Constructor" begin
        usm6 = USM6{Float64}(1.0, 0.1, 0.2, 0.3, 0.4, 0.5)
        @test usm6 isa USM6{Float64}
        @test usm6.C == 1.0
        @test usm6.Rf1 == 0.1
        @test usm6.Rf2 == 0.2
        @test usm6.σ1 == 0.3
        @test usm6.σ2 == 0.4
        @test usm6.σ3 == 0.5
    end
    
    @testset "USM6 - Constructor from Vector" begin
        vec = [1.0, 0.1, 0.2, 0.3, 0.4, 0.5]
        usm6 = USM6(vec)
        @test usm6 isa USM6{Float64}
        @test usm6.C == 1.0
        @test usm6.σ3 == 0.5
    end
    
    @testset "USM6 - Type Promotion" begin
        usm6 = USM6(1.0, 0.1f0, 0.2, 0, 0.4, 0.5)
        @test usm6 isa USM6
        @test eltype(usm6) <: AbstractFloat
    end
    
    @testset "USM6 - StaticVector Constructor" begin
        svec = SVector{6}(1.0, 0.1, 0.2, 0.3, 0.4, 0.5)
        usm6 = USM6(svec)
        @test usm6 isa USM6{Float64}
    end
    
    @testset "USM6 - params() Method" begin
        usm6 = USM6(1.0, 0.1, 0.2, 0.3, 0.4, 0.5)
        p = params(usm6)
        
        @test p isa SVector{6, Float64}
        @test p[1] == 1.0
        @test p[6] == 0.5
    end
    
    @testset "USM6 - Base.one()" begin
        usm6_one = one(USM6)
        @test usm6_one isa USM6{Float64}
        @test all(params(usm6_one) .== 0.0)
    end
    
    @testset "USM6 - Base.getindex" begin
        usm6 = USM6(1.0, 0.1, 0.2, 0.3, 0.4, 0.5)
        
        @test usm6[1] == 1.0
        @test usm6[2] == 0.1
        @test usm6[3] == 0.2
        @test usm6[4] == 0.3
        @test usm6[5] == 0.4
        @test usm6[6] == 0.5
        
        @test_throws BoundsError usm6[0]
        @test_throws BoundsError usm6[7]
    end
    
    @testset "USMEM - Inner Constructor" begin
        usmem = USMEM{Float64}(1.0, 0.1, 0.2, 0.3, 0.4, 0.5)
        @test usmem isa USMEM{Float64}
        @test usmem.C == 1.0
        @test usmem.Rf1 == 0.1
        @test usmem.Rf2 == 0.2
        @test usmem.a1 == 0.3
        @test usmem.a2 == 0.4
        @test usmem.a3 == 0.5
    end
    
    @testset "USMEM - Constructor from Vector" begin
        vec = [1.0, 0.1, 0.2, 0.3, 0.4, 0.5]
        usmem = USMEM(vec)
        @test usmem isa USMEM{Float64}
        @test usmem.a3 == 0.5
    end
    
    @testset "USMEM - Type Promotion" begin
        usmem = USMEM(1.0, 0.1f0, 0.2, 0, 0.4, 0.5)
        @test usmem isa USMEM
        @test eltype(usmem) <: AbstractFloat
    end
    
    @testset "USMEM - params() Method" begin
        usmem = USMEM(1.0, 0.1, 0.2, 0.3, 0.4, 0.5)
        p = params(usmem)
        
        @test p isa SVector{6, Float64}
        @test p[1] == 1.0
        @test p[6] == 0.5
    end
    
    @testset "USMEM - Base.one()" begin
        usmem_one = one(USMEM)
        @test usmem_one isa USMEM{Float64}
        @test all(params(usmem_one) .== 0.0)
    end
    
    @testset "USMEM - Base.getindex" begin
        usmem = USMEM(1.0, 0.1, 0.2, 0.3, 0.4, 0.5)
        
        @test usmem[1] == 1.0
        @test usmem[2] == 0.1
        @test usmem[3] == 0.2
        @test usmem[4] == 0.3
        @test usmem[5] == 0.4
        @test usmem[6] == 0.5
        
        @test_throws BoundsError usmem[0]
        @test_throws BoundsError usmem[7]
    end
    
    @testset "USM7 ↔ Keplerian Transformations" begin
        μ = 398600.4418
        
        # Create Keplerian elements (circular orbit)
        kep = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.5)
        
        # Convert to USM7
        usm7 = USM7(kep, μ)
        @test usm7 isa USM7
        @test all(isfinite.(params(usm7)))
        
        # Convert back to Keplerian
        kep_back = Keplerian(usm7, μ)
        @test kep_back isa Keplerian
        @test all(isfinite.(params(kep_back)))
        
        # Check round-trip accuracy
        @test kep_back.a ≈ kep.a atol=1e-6 rtol=1e-8
        @test kep_back.e ≈ kep.e atol=1e-8 rtol=1e-8
        @test kep_back.i ≈ kep.i atol=1e-8 rtol=1e-8
    end
    
    @testset "USM7 ↔ USM6 Transformations" begin
        μ = 398600.4418
        
        # Create USM7
        usm7 = USM7(1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
        
        # Convert to USM6
        usm6 = USM6(usm7, μ)
        @test usm6 isa USM6
        @test all(isfinite.(params(usm6)))
        
        # Convert back to USM7
        usm7_back = USM7(usm6, μ)
        @test usm7_back isa USM7
        @test all(isfinite.(params(usm7_back)))
        
        # Check round-trip accuracy
        @test usm7_back.C ≈ usm7.C atol=1e-8 rtol=1e-8
        @test usm7_back.Rf1 ≈ usm7.Rf1 atol=1e-8 rtol=1e-8
        @test usm7_back.Rf2 ≈ usm7.Rf2 atol=1e-8 rtol=1e-8
    end
    
    @testset "USM7 ↔ USMEM Transformations" begin
        μ = 398600.4418
        
        # Create USM7
        usm7 = USM7(1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
        
        # Convert to USMEM
        usmem = USMEM(usm7, μ)
        @test usmem isa USMEM
        @test all(isfinite.(params(usmem)))
        
        # Convert back to USM7
        usm7_back = USM7(usmem, μ)
        @test usm7_back isa USM7
        @test all(isfinite.(params(usm7_back)))
    end
    
    @testset "Singularity-Free Properties" begin
        μ = 398600.4418
        
        # Test circular orbit (e=0)
        kep_circ = Keplerian(7000.0, 0.0, 0.1, 0.0, 0.0, 0.5)
        usm7_circ = USM7(kep_circ, μ)
        @test all(isfinite.(params(usm7_circ)))
        @test all(.!isnan.(params(usm7_circ)))
        
        # Test equatorial orbit (i=0)
        kep_equat = Keplerian(7000.0, 0.1, 0.0, 0.0, 0.0, 0.5)
        usm7_equat = USM7(kep_equat, μ)
        @test all(isfinite.(params(usm7_equat)))
        @test all(.!isnan.(params(usm7_equat)))
        
        # Test circular equatorial (e=0, i=0)
        kep_ce = Keplerian(7000.0, 0.0, 0.0, 0.0, 0.0, 0.5)
        usm7_ce = USM7(kep_ce, μ)
        @test all(isfinite.(params(usm7_ce)))
        @test all(.!isnan.(params(usm7_ce)))
    end
    
    @testset "Numerical Stability - Various Orbit Types" begin
        μ = 398600.4418
        
        # Highly elliptical orbit
        kep_ellip = Keplerian(10000.0, 0.7, 0.5, 1.0, 1.5, 2.0)
        usm7_ellip = USM7(kep_ellip, μ)
        @test all(isfinite.(params(usm7_ellip)))
        kep_ellip_back = Keplerian(usm7_ellip, μ)
        @test kep_ellip_back.a ≈ kep_ellip.a atol=1e-6 rtol=1e-8
        
        # High inclination orbit
        kep_incl = Keplerian(7000.0, 0.1, 1.5, 0.5, 1.0, 0.5)
        usm7_incl = USM7(kep_incl, μ)
        @test all(isfinite.(params(usm7_incl)))
        
        # Retrograde orbit
        kep_retro = Keplerian(7000.0, 0.1, 2.5, 0.5, 1.0, 0.5)
        usm7_retro = USM7(kep_retro, μ)
        @test all(isfinite.(params(usm7_retro)))
    end
    
    @testset "Type Stability" begin
        μ = 398600.4418
        
        # Float64
        kep_f64 = Keplerian(7000.0, 0.1, 0.1, 0.2, 0.3, 0.5)
        usm7_f64 = USM7(kep_f64, μ)
        @test eltype(usm7_f64) === Float64
        
        # Float32
        kep_f32 = Keplerian(7000.0f0, 0.1f0, 0.1f0, 0.2f0, 0.3f0, 0.5f0)
        usm7_f32 = USM7(kep_f32, Float32(μ))
        @test eltype(usm7_f32) <: AbstractFloat
    end
end
