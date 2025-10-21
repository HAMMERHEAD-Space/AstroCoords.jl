using Test
using AstroCoords
using StaticArrays

@testset "Spherical Constructor Coverage" begin
    # Target: Lines 26, 35-37, 40-41, 43, 46-47
    # These are constructor variants - no literature needed (API tests)
    
    @testset "AbstractVector constructor" begin
        # Line 26: Spherical(X::AbstractVector{T})
        vec = [7000.0, 0.5, 0.3, 1.5, 0.01, 0.02]
        sph = Spherical(vec)
        @test sph isa Spherical{Float64}
        @test sph.r ≈ 7000.0
        @test sph.θ ≈ 0.5
        @test sph.ϕ ≈ 0.3
        @test sph.ṙ ≈ 1.5
        @test sph.θdot ≈ 0.01
        @test sph.ϕdot ≈ 0.02
    end
    
    @testset "Individual args with type promotion" begin
        # Lines 35-37: Spherical(r::R, θ::T, ϕ::P, ṙ::RD, θdot::TD, ϕdot::PD) with promote_type
        # Test Int → Float64 promotion
        sph1 = Spherical(7000, 0.5, 0.3, 1.5, 0.01, 0.02)
        @test sph1 isa Spherical{Float64}
        @test sph1.r ≈ 7000.0
        
        # Test Float32 → Float64 promotion
        sph2 = Spherical(7000.0f0, 0.5f0, 0.3, 1.5, 0.01, 0.02)
        @test sph2 isa Spherical{Float64}
        @test sph2.r ≈ 7000.0
        
        # Test mixed Int/Float32/Float64
        sph3 = Spherical(7000, 0.5f0, 0.3, 1.5f0, 0.01, 0.02f0)
        @test sph3 isa Spherical{Float64}
        @test sph3.θ ≈ 0.5
    end
    
    @testset "StaticVector constructor" begin
        # Lines 40-41, 43: Spherical(g::StaticVector{N,T})
        sv = SVector{6,Float64}(7000.0, 0.5, 0.3, 1.5, 0.01, 0.02)
        sph = Spherical(sv)
        @test sph isa Spherical{Float64}
        @test sph.r ≈ 7000.0
        @test sph.ϕdot ≈ 0.02
        
        # Test with Float32 StaticVector
        sv32 = SVector{6,Float32}(7000.0f0, 0.5f0, 0.3f0, 1.5f0, 0.01f0, 0.02f0)
        sph32 = Spherical(sv32)
        @test sph32 isa Spherical{Float32}
        @test sph32.r ≈ 7000.0f0
    end
    
    @testset "Spherical{T}(StaticVector) type conversion" begin
        # Line 43: Spherical{T}(g::StaticVector)
        sv = SVector{6,Float32}(7000.0f0, 0.5f0, 0.3f0, 1.5f0, 0.01f0, 0.02f0)
        sph = Spherical{Float64}(sv)
        @test sph isa Spherical{Float64}
        @test sph.r ≈ 7000.0
    end
    
    @testset "Base.one() with default type" begin
        # Lines 46-47: Base.one(::Type{S}; T::DataType=Float64)
        sph_one = Base.one(Spherical)
        @test sph_one isa Spherical{Float64}
        @test sph_one.r == 0.0
        @test sph_one.θ == 0.0
        @test sph_one.ϕ == 0.0
        @test sph_one.ṙ == 0.0
        @test sph_one.θdot == 0.0
        @test sph_one.ϕdot == 0.0
    end
    
    @testset "Base.one() with custom type" begin
        # Lines 46-47: Base.one with T=Float32
        sph_one32 = Base.one(Spherical; T=Float32)
        @test sph_one32 isa Spherical{Float32}
        @test sph_one32.r == 0.0f0
        @test sph_one32.ṙ == 0.0f0
    end
    
    @testset "Field access via indexing" begin
        # Lines 50-69: Base.getindex(p::Spherical{T}, i::Int)
        sph = Spherical(7000.0, 0.5, 0.3, 1.5, 0.01, 0.02)
        @test sph[1] ≈ 7000.0  # r
        @test sph[2] ≈ 0.5     # θ
        @test sph[3] ≈ 0.3     # ϕ
        @test sph[4] ≈ 1.5     # ṙ
        @test sph[5] ≈ 0.01    # θdot
        @test sph[6] ≈ 0.02    # ϕdot
        @test_throws BoundsError sph[7]
        @test_throws BoundsError sph[0]
    end
    
    @testset "params() extraction" begin
        # Line 45: params(g::Spherical{T})
        sph = Spherical(7000.0, 0.5, 0.3, 1.5, 0.01, 0.02)
        p = params(sph)
        @test p isa SVector{6,Float64}
        @test p[1] ≈ 7000.0
        @test p[6] ≈ 0.02
    end
    
    @testset "Boundary cases" begin
        # Test edge cases: r=0, θ=0/π, ϕ=0/2π, zero velocities
        sph_zero_r = Spherical(0.0, 0.5, 0.3, 0.0, 0.0, 0.0)
        @test sph_zero_r.r == 0.0
        @test sph_zero_r.ṙ == 0.0
        
        sph_zero_theta = Spherical(7000.0, 0.0, 0.3, 1.5, 0.0, 0.02)
        @test sph_zero_theta.θ == 0.0
        @test sph_zero_theta.θdot == 0.0
        
        sph_pi_theta = Spherical(7000.0, π, 0.3, 1.5, 0.01, 0.02)
        @test sph_pi_theta.θ ≈ π
        
        sph_zero_phi = Spherical(7000.0, 0.5, 0.0, 1.5, 0.01, 0.0)
        @test sph_zero_phi.ϕ == 0.0
        @test sph_zero_phi.ϕdot == 0.0
        
        sph_twopi_phi = Spherical(7000.0, 0.5, 2π, 1.5, 0.01, 0.02)
        @test sph_twopi_phi.ϕ ≈ 2π
        
        sph_zero_vel = Spherical(7000.0, 0.5, 0.3, 0.0, 0.0, 0.0)
        @test sph_zero_vel.ṙ == 0.0
        @test sph_zero_vel.θdot == 0.0
        @test sph_zero_vel.ϕdot == 0.0
    end
end
