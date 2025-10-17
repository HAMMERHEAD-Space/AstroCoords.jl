using Test
using AstroCoords
using StaticArrays

@testset "Keplerian Coordinate Set" begin
    @testset "Inner constructor Keplerian{T}" begin
        k = Keplerian{Float64}(7000.0, 0.01, 0.1, 0.2, 0.3, 0.4)
        @test k isa Keplerian{Float64}
        @test k.a ≈ 7000.0
        @test k.e ≈ 0.01
        @test k.i ≈ 0.1
        @test k.Ω ≈ 0.2
        @test k.ω ≈ 0.3
        @test k.f ≈ 0.4
    end

    @testset "Constructor with type promotion" begin
        # Float32 + Float64 should promote to Float64
        k = Keplerian(7000.0f0, 0.01, 0.1f0, 0.2, 0.3f0, 0.4)
        @test k isa Keplerian{Float64}
        @test k.a ≈ 7000.0
        @test k.e ≈ 0.01
        @test k.i ≈ 0.1
        @test k.Ω ≈ 0.2
        @test k.ω ≈ 0.3
        @test k.f ≈ 0.4
    end

    @testset "StaticVector constructor" begin
        vec = SVector{6}(7000.0, 0.01, 0.1, 0.2, 0.3, 0.4)
        k = Keplerian(vec)
        @test k isa Keplerian{Float64}
        @test k.a ≈ 7000.0
        @test k.e ≈ 0.01
        @test k.i ≈ 0.1
        @test k.Ω ≈ 0.2
        @test k.ω ≈ 0.3
        @test k.f ≈ 0.4
    end

    @testset "Constructor from vector" begin
        vec = [7000.0, 0.01, 0.1, 0.2, 0.3, 0.4]
        k = Keplerian(vec)
        @test k isa Keplerian{Float64}
        @test k.a ≈ 7000.0
        @test k.e ≈ 0.01
        @test k.i ≈ 0.1
        @test k.Ω ≈ 0.2
        @test k.ω ≈ 0.3
        @test k.f ≈ 0.4
    end

    @testset "Base.one(Keplerian)" begin
        k_one = Base.one(Keplerian)
        @test k_one isa Keplerian{Float64}
        @test k_one.a ≈ 0.0
        @test k_one.e ≈ 0.0
        @test k_one.i ≈ 0.0
        @test k_one.Ω ≈ 0.0
        @test k_one.ω ≈ 0.0
        @test k_one.f ≈ 0.0

        k_one_f32 = Base.one(Keplerian; T=Float32)
        @test k_one_f32 isa Keplerian{Float32}
        @test k_one_f32.a ≈ 0.0f0
    end

    @testset "Base.getindex for all 6 indices" begin
        k = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.4)
        
        @test k[1] ≈ 7000.0  # a
        @test k[2] ≈ 0.01    # e
        @test k[3] ≈ 0.1     # i
        @test k[4] ≈ 0.2     # Ω
        @test k[5] ≈ 0.3     # ω
        @test k[6] ≈ 0.4     # f
        
        @test_throws BoundsError k[0]
        @test_throws BoundsError k[7]
    end

    @testset "Various element ranges" begin
        # Test with different semi-major axes
        k_low = Keplerian(6500.0, 0.01, 0.1, 0.2, 0.3, 0.4)
        k_high = Keplerian(42164.0, 0.01, 0.1, 0.2, 0.3, 0.4)  # GEO altitude
        @test k_low.a < k_high.a
        
        # Test with different eccentricities
        k_circ = Keplerian(7000.0, 0.0, 0.1, 0.2, 0.3, 0.4)  # Circular
        k_ell = Keplerian(7000.0, 0.5, 0.1, 0.2, 0.3, 0.4)   # Elliptical
        @test k_circ.e < k_ell.e
        
        # Test with different inclinations
        k_eq = Keplerian(7000.0, 0.01, 0.0, 0.2, 0.3, 0.4)      # Equatorial
        k_polar = Keplerian(7000.0, 0.01, π/2, 0.2, 0.3, 0.4)   # Polar
        @test k_eq.i < k_polar.i
    end

    @testset "Singularity cases" begin
        # Circular orbit (e=0)
        k_circ = Keplerian(7000.0, 0.0, 0.5, 1.0, 0.0, 0.0)
        @test k_circ.e ≈ 0.0
        @test k_circ.ω ≈ 0.0  # ω undefined for circular, typically set to 0
        
        # Equatorial orbit (i=0)
        k_eq = Keplerian(7000.0, 0.1, 0.0, 0.0, 0.5, 0.3)
        @test k_eq.i ≈ 0.0
        @test k_eq.Ω ≈ 0.0  # Ω undefined for equatorial, typically set to 0
        
        # Circular equatorial orbit (e=0, i=0)
        k_circ_eq = Keplerian(7000.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        @test k_circ_eq.e ≈ 0.0
        @test k_circ_eq.i ≈ 0.0
        @test k_circ_eq.Ω ≈ 0.0
        @test k_circ_eq.ω ≈ 0.0
    end

    @testset "Edge cases: e=1, e>1" begin
        # Parabolic orbit (e=1)
        k_para = Keplerian(7000.0, 1.0, 0.5, 1.0, 0.5, 0.3)
        @test k_para.e ≈ 1.0
        @test k_para.a > 0.0  # Note: For e=1, semi-major axis is infinite, but we store finite value
        
        # Hyperbolic orbit (e>1)
        k_hyper = Keplerian(-10000.0, 1.5, 0.5, 1.0, 0.5, 0.3)  # Negative a for hyperbola
        @test k_hyper.e > 1.0
        @test k_hyper.a < 0.0  # Negative semi-major axis for hyperbolic orbit
        
        # Near-parabolic (e slightly less than 1)
        k_near_para = Keplerian(7000.0, 0.9999, 0.5, 1.0, 0.5, 0.3)
        @test k_near_para.e < 1.0
        @test k_near_para.e > 0.999
    end

    @testset "True anomaly full range" begin
        # Test various true anomaly values (0 to 2π)
        for f in [0.0, π/4, π/2, π, 3π/2, 2π - 0.01]
            k = Keplerian(7000.0, 0.1, 0.5, 1.0, 0.5, f)
            @test k.f ≈ f rtol=1e-10
        end
    end

    @testset "Type stability" begin
        k_f64 = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.4)
        @test typeof(k_f64.a) === Float64
        @test typeof(k_f64.e) === Float64
        
        k_f32 = Keplerian(7000.0f0, 0.01f0, 0.1f0, 0.2f0, 0.3f0, 0.4f0)
        @test typeof(k_f32.a) === Float32
        @test typeof(k_f32.e) === Float32
    end
end
