using Test
using AstroCoords
using StaticArrays

@testset "Spherical and Cylindrical Coordinate Sets" begin
    @testset "Spherical - Inner Constructor" begin
        sph = Spherical{Float64}(7000.0, 0.5, 1.0, 1.0, 0.001, 0.002)
        @test sph isa Spherical{Float64}
        @test sph.r == 7000.0
        @test sph.θ == 0.5
        @test sph.ϕ == 1.0
        @test sph.ṙ == 1.0
        @test sph.θdot == 0.001
        @test sph.ϕdot == 0.002
    end
    
    @testset "Spherical - Constructor from Vector" begin
        vec = [7000.0, 0.5, 1.0, 1.0, 0.001, 0.002]
        sph = Spherical(vec)
        @test sph isa Spherical{Float64}
        @test sph.r == 7000.0
        @test sph.ϕdot == 0.002
        
        # Test with Int vector
        vec_int = [7000, 0, 1, 1, 0, 0]
        sph_int = Spherical(vec_int)
        @test sph_int isa Spherical{Int}
    end
    
    @testset "Spherical - Type Promotion" begin
        # Mixed types
        sph = Spherical(7000.0, 0.5f0, 1.0, 1, 0.001, 0.002)
        @test sph isa Spherical
        @test eltype(sph) <: AbstractFloat
        
        # All Float32
        sph_f32 = Spherical(7000.0f0, 0.5f0, 1.0f0, 1.0f0, 0.001f0, 0.002f0)
        @test sph_f32 isa Spherical{Float32}
    end
    
    @testset "Spherical - StaticVector Constructor" begin
        svec = SVector{6}(7000.0, 0.5, 1.0, 1.0, 0.001, 0.002)
        sph = Spherical(svec)
        @test sph isa Spherical{Float64}
        @test sph.r == 7000.0
    end
    
    @testset "Spherical - params() Method" begin
        sph = Spherical(7000.0, 0.5, 1.0, 1.0, 0.001, 0.002)
        p = params(sph)
        
        @test p isa SVector{6, Float64}
        @test p[1] == 7000.0
        @test p[2] == 0.5
        @test p[3] == 1.0
        @test p[4] == 1.0
        @test p[5] == 0.001
        @test p[6] == 0.002
    end
    
    @testset "Spherical - Base.one()" begin
        sph_one = one(Spherical)
        @test sph_one isa Spherical{Float64}
        @test all(params(sph_one) .== 0.0)
        
        sph_one_f32 = one(Spherical, T=Float32)
        @test sph_one_f32 isa Spherical{Float32}
    end
    
    @testset "Spherical - Base.getindex" begin
        sph = Spherical(7000.0, 0.5, 1.0, 1.0, 0.001, 0.002)
        
        @test sph[1] == 7000.0  # r
        @test sph[2] == 0.5     # θ
        @test sph[3] == 1.0     # ϕ
        @test sph[4] == 1.0     # ṙ
        @test sph[5] == 0.001   # θdot
        @test sph[6] == 0.002   # ϕdot
        
        @test_throws BoundsError sph[0]
        @test_throws BoundsError sph[7]
    end
    
    @testset "Spherical ↔ Cartesian Transformations" begin
        μ = 398600.4418
        
        # Create Spherical coordinates
        sph = Spherical(7000.0, π/4, π/6, 1.0, 0.001, 0.001)
        
        # Convert to Cartesian
        cart = Cartesian(sph, μ)
        @test cart isa Cartesian
        @test all(isfinite.(params(cart)))
        
        # Convert back to Spherical
        sph_back = Spherical(cart, μ)
        @test sph_back isa Spherical
        @test all(isfinite.(params(sph_back)))
    end
    
    @testset "Spherical Round-Trip" begin
        μ = 398600.4418
        
        sph_orig = Spherical(7000.0, π/4, π/6, 1.0, 0.001, 0.001)
        cart = Cartesian(sph_orig, μ)
        sph_back = Spherical(cart, μ)
        
        @test sph_back.r ≈ sph_orig.r atol=1e-6 rtol=1e-8
        @test sph_back.θ ≈ sph_orig.θ atol=1e-8 rtol=1e-8
        @test sph_back.ϕ ≈ sph_orig.ϕ atol=1e-8 rtol=1e-8
    end
    
    @testset "Spherical - Coordinate Ranges" begin
        μ = 398600.4418
        
        # Test r > 0
        sph_large_r = Spherical(20000.0, π/4, π/3, 1.0, 0.001, 0.001)
        cart_large = Cartesian(sph_large_r, μ)
        @test all(isfinite.(params(cart_large)))
        
        # Test θ = 0 (north pole)
        sph_pole_n = Spherical(7000.0, 0.0, 0.0, 1.0, 0.001, 0.001)
        cart_pole_n = Cartesian(sph_pole_n, μ)
        @test all(isfinite.(params(cart_pole_n)))
        
        # Test θ = π (south pole)
        sph_pole_s = Spherical(7000.0, π, 0.0, 1.0, 0.001, 0.001)
        cart_pole_s = Cartesian(sph_pole_s, μ)
        @test all(isfinite.(params(cart_pole_s)))
        
        # Test ϕ = 0
        sph_phi0 = Spherical(7000.0, π/4, 0.0, 1.0, 0.001, 0.001)
        cart_phi0 = Cartesian(sph_phi0, μ)
        @test all(isfinite.(params(cart_phi0)))
        
        # Test ϕ near 2π
        sph_phi2pi = Spherical(7000.0, π/4, 2π - 0.01, 1.0, 0.001, 0.001)
        cart_phi2pi = Cartesian(sph_phi2pi, μ)
        @test all(isfinite.(params(cart_phi2pi)))
    end
    
    @testset "Cylindrical - Inner Constructor" begin
        cyl = Cylindrical{Float64}(5000.0, 0.5, 2000.0, 1.0, 0.001, 0.5)
        @test cyl isa Cylindrical{Float64}
        @test cyl.ρ == 5000.0
        @test cyl.θ == 0.5
        @test cyl.z == 2000.0
        @test cyl.ρdot == 1.0
        @test cyl.θdot == 0.001
        @test cyl.ż == 0.5
    end
    
    @testset "Cylindrical - Constructor from Vector" begin
        vec = [5000.0, 0.5, 2000.0, 1.0, 0.001, 0.5]
        cyl = Cylindrical(vec)
        @test cyl isa Cylindrical{Float64}
        @test cyl.ρ == 5000.0
        @test cyl.ż == 0.5
    end
    
    @testset "Cylindrical - Type Promotion" begin
        cyl = Cylindrical(5000.0, 0.5f0, 2000.0, 1, 0.001, 0.5)
        @test cyl isa Cylindrical
        @test eltype(cyl) <: AbstractFloat
    end
    
    @testset "Cylindrical - params() Method" begin
        cyl = Cylindrical(5000.0, 0.5, 2000.0, 1.0, 0.001, 0.5)
        p = params(cyl)
        
        @test p isa SVector{6, Float64}
        @test p[1] == 5000.0
        @test p[6] == 0.5
    end
    
    @testset "Cylindrical - Base.one()" begin
        cyl_one = one(Cylindrical)
        @test cyl_one isa Cylindrical{Float64}
        @test all(params(cyl_one) .== 0.0)
    end
    
    @testset "Cylindrical - Base.getindex" begin
        cyl = Cylindrical(5000.0, 0.5, 2000.0, 1.0, 0.001, 0.5)
        
        @test cyl[1] == 5000.0  # ρ
        @test cyl[2] == 0.5     # θ
        @test cyl[3] == 2000.0  # z
        @test cyl[4] == 1.0     # ρdot
        @test cyl[5] == 0.001   # θdot
        @test cyl[6] == 0.5     # ż
        
        @test_throws BoundsError cyl[0]
        @test_throws BoundsError cyl[7]
    end
    
    @testset "Cylindrical ↔ Cartesian Transformations" begin
        μ = 398600.4418
        
        cyl = Cylindrical(5000.0, π/4, 2000.0, 1.0, 0.001, 0.5)
        
        cart = Cartesian(cyl, μ)
        @test cart isa Cartesian
        @test all(isfinite.(params(cart)))
        
        cyl_back = Cylindrical(cart, μ)
        @test cyl_back isa Cylindrical
        @test all(isfinite.(params(cyl_back)))
    end
    
    @testset "Cylindrical Round-Trip" begin
        μ = 398600.4418
        
        cyl_orig = Cylindrical(5000.0, π/4, 2000.0, 1.0, 0.001, 0.5)
        cart = Cartesian(cyl_orig, μ)
        cyl_back = Cylindrical(cart, μ)
        
        @test cyl_back.ρ ≈ cyl_orig.ρ atol=1e-6 rtol=1e-8
        @test cyl_back.θ ≈ cyl_orig.θ atol=1e-8 rtol=1e-8
        @test cyl_back.z ≈ cyl_orig.z atol=1e-6 rtol=1e-8
    end
    
    @testset "Type Stability" begin
        μ = 398600.4418
        
        # Spherical Float64
        sph_f64 = Spherical(7000.0, π/4, π/6, 1.0, 0.001, 0.001)
        cart_f64 = Cartesian(sph_f64, μ)
        @test eltype(cart_f64) === Float64
        
        # Spherical Float32
        sph_f32 = Spherical(7000.0f0, Float32(π/4), Float32(π/6), 1.0f0, 0.001f0, 0.001f0)
        cart_f32 = Cartesian(sph_f32, Float32(μ))
        @test eltype(cart_f32) <: AbstractFloat
        
        # Cylindrical Float64
        cyl_f64 = Cylindrical(5000.0, π/4, 2000.0, 1.0, 0.001, 0.5)
        cart_cyl_f64 = Cartesian(cyl_f64, μ)
        @test eltype(cart_cyl_f64) === Float64
    end
end
