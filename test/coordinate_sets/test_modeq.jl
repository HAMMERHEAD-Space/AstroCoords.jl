using Test
using AstroCoords
using StaticArrays

@testset "ModEq and J2EqOE Coordinate Sets" begin
    @testset "ModEq - Inner Constructor" begin
        modeq = ModEq{Float64}(7000.0, 0.1, 0.2, 0.05, 0.06, 1.0)
        @test modeq isa ModEq{Float64}
        @test modeq.p == 7000.0
        @test modeq.f == 0.1
        @test modeq.g == 0.2
        @test modeq.h == 0.05
        @test modeq.k == 0.06
        @test modeq.L == 1.0
    end
    
    @testset "ModEq - Constructor from Vector" begin
        vec = [7000.0, 0.1, 0.2, 0.05, 0.06, 1.0]
        modeq = ModEq(vec)
        @test modeq isa ModEq{Float64}
        @test modeq.p == 7000.0
        @test modeq.L == 1.0
        
        # Test with Int vector
        vec_int = [7000, 0, 0, 0, 0, 1]
        modeq_int = ModEq(vec_int)
        @test modeq_int isa ModEq{Int}
    end
    
    @testset "ModEq - Type Promotion" begin
        # Mixed types
        modeq = ModEq(7000.0, 0.1f0, 0.2, 0, 0.06, 1.0)
        @test modeq isa ModEq
        @test eltype(modeq) <: AbstractFloat
        
        # All Float32
        modeq_f32 = ModEq(7000.0f0, 0.1f0, 0.2f0, 0.05f0, 0.06f0, 1.0f0)
        @test modeq_f32 isa ModEq{Float32}
    end
    
    @testset "ModEq - StaticVector Constructor" begin
        svec = SVector{6}(7000.0, 0.1, 0.2, 0.05, 0.06, 1.0)
        modeq = ModEq(svec)
        @test modeq isa ModEq{Float64}
        @test modeq.p == 7000.0
    end
    
    @testset "ModEq - Constructor from Vector (line 38)" begin
        # Test the specific constructor that takes AbstractVector
        vec = [7000.0, 0.1, 0.2, 0.05, 0.06, 1.0]
        modeq = ModEq(vec)
        @test modeq isa ModEq{Float64}
        @test length(params(modeq)) == 6
    end
    
    @testset "ModEq - Base.one()" begin
        modeq_one = one(ModEq)
        @test modeq_one isa ModEq{Float64}
        @test all(params(modeq_one) .== 0.0)
        
        modeq_one_f32 = one(ModEq, T=Float32)
        @test modeq_one_f32 isa ModEq{Float32}
    end
    
    @testset "ModEq - params() Method" begin
        modeq = ModEq(7000.0, 0.1, 0.2, 0.05, 0.06, 1.0)
        p = params(modeq)
        
        @test p isa SVector{6, Float64}
        @test p[1] == 7000.0
        @test p[2] == 0.1
        @test p[3] == 0.2
        @test p[4] == 0.05
        @test p[5] == 0.06
        @test p[6] == 1.0
    end
    
    @testset "ModEq - Base.getindex" begin
        modeq = ModEq(7000.0, 0.1, 0.2, 0.05, 0.06, 1.0)
        
        @test modeq[1] == 7000.0  # p
        @test modeq[2] == 0.1     # f
        @test modeq[3] == 0.2     # g
        @test modeq[4] == 0.05    # h
        @test modeq[5] == 0.06    # k
        @test modeq[6] == 1.0     # L
        
        @test_throws BoundsError modeq[0]
        @test_throws BoundsError modeq[7]
    end
    
    @testset "ModEq ↔ Keplerian Transformations" begin
        μ = 398600.4418
        
        # Create Keplerian elements
        kep = Keplerian(7000.0, 0.1, 0.2, 0.3, 0.4, 0.5)
        
        # Convert to ModEq
        modeq = ModEq(kep, μ)
        @test modeq isa ModEq
        @test all(isfinite.(params(modeq)))
        
        # Convert back to Keplerian
        kep_back = Keplerian(modeq, μ)
        @test kep_back isa Keplerian
        @test all(isfinite.(params(kep_back)))
    end
    
    @testset "ModEq Round-Trip" begin
        μ = 398600.4418
        
        kep_orig = Keplerian(7000.0, 0.1, 0.2, 0.3, 0.4, 0.5)
        modeq = ModEq(kep_orig, μ)
        kep_back = Keplerian(modeq, μ)
        
        @test kep_back.a ≈ kep_orig.a atol=1e-6 rtol=1e-8
        @test kep_back.e ≈ kep_orig.e atol=1e-8 rtol=1e-8
        @test kep_back.i ≈ kep_orig.i atol=1e-8 rtol=1e-8
        @test kep_back.Ω ≈ kep_orig.Ω atol=1e-8 rtol=1e-8
        @test kep_back.ω ≈ kep_orig.ω atol=1e-8 rtol=1e-8
    end
    
    @testset "ModEq - Different Orbit Types" begin
        μ = 398600.4418
        
        # Circular orbit (e ≈ 0)
        kep_circ = Keplerian(7000.0, 0.001, 0.2, 0.3, 0.4, 0.5)
        modeq_circ = ModEq(kep_circ, μ)
        @test all(isfinite.(params(modeq_circ)))
        @test all(.!isnan.(params(modeq_circ)))
        
        # Elliptical orbit
        kep_ellip = Keplerian(10000.0, 0.3, 0.5, 1.0, 1.5, 2.0)
        modeq_ellip = ModEq(kep_ellip, μ)
        @test all(isfinite.(params(modeq_ellip)))
        
        # Near-parabolic orbit (e ≈ 1)
        kep_para = Keplerian(7000.0, 0.95, 0.2, 0.3, 0.4, 0.5)
        modeq_para = ModEq(kep_para, μ)
        @test all(isfinite.(params(modeq_para)))
    end
    
    @testset "J2EqOE - Inner Constructor" begin
        j2eq = J2EqOE{Float64}(0.001, 0.1, 0.2, 0.05, 0.06, 1.0)
        @test j2eq isa J2EqOE{Float64}
        @test j2eq.n == 0.001
        @test j2eq.h == 0.1
        @test j2eq.k == 0.2
        @test j2eq.p == 0.05
        @test j2eq.q == 0.06
        @test j2eq.L == 1.0
    end
    
    @testset "J2EqOE - Constructor from Vector" begin
        vec = [0.001, 0.1, 0.2, 0.05, 0.06, 1.0]
        j2eq = J2EqOE(vec)
        @test j2eq isa J2EqOE{Float64}
        @test j2eq.n == 0.001
        @test j2eq.L == 1.0
    end
    
    @testset "J2EqOE - Type Promotion" begin
        j2eq = J2EqOE(0.001, 0.1f0, 0.2, 0, 0.06, 1.0)
        @test j2eq isa J2EqOE
        @test eltype(j2eq) <: AbstractFloat
    end
    
    @testset "J2EqOE - params() Method" begin
        j2eq = J2EqOE(0.001, 0.1, 0.2, 0.05, 0.06, 1.0)
        p = params(j2eq)
        
        @test p isa SVector{6, Float64}
        @test p[1] == 0.001
        @test p[6] == 1.0
    end
    
    @testset "J2EqOE - Base.one()" begin
        j2eq_one = one(J2EqOE)
        @test j2eq_one isa J2EqOE{Float64}
        @test all(params(j2eq_one) .== 0.0)
    end
    
    @testset "J2EqOE - Base.getindex" begin
        j2eq = J2EqOE(0.001, 0.1, 0.2, 0.05, 0.06, 1.0)
        
        @test j2eq[1] == 0.001  # n
        @test j2eq[2] == 0.1    # h
        @test j2eq[3] == 0.2    # k
        @test j2eq[4] == 0.05   # p
        @test j2eq[5] == 0.06   # q
        @test j2eq[6] == 1.0    # L
        
        @test_throws BoundsError j2eq[0]
        @test_throws BoundsError j2eq[7]
    end
    
    @testset "J2EqOE ↔ Cartesian Transformations" begin
        μ = 398600.4418
        J2 = 1.08263e-3
        Req = 6378.137
        
        # Create J2EqOE
        j2eq = J2EqOE(0.001, 0.01, 0.02, 0.005, 0.006, 1.0)
        
        # Convert to Cartesian
        cart = Cartesian(j2eq, μ, J2, Req)
        @test cart isa Cartesian
        @test all(isfinite.(params(cart)))
        
        # Convert back to J2EqOE
        j2eq_back = J2EqOE(cart, μ, J2, Req)
        @test j2eq_back isa J2EqOE
        @test all(isfinite.(params(j2eq_back)))
    end
    
    @testset "Type Stability" begin
        μ = 398600.4418
        
        # ModEq Float64
        kep_f64 = Keplerian(7000.0, 0.1, 0.2, 0.3, 0.4, 0.5)
        modeq_f64 = ModEq(kep_f64, μ)
        @test eltype(modeq_f64) === Float64
        
        # ModEq Float32
        kep_f32 = Keplerian(7000.0f0, 0.1f0, 0.2f0, 0.3f0, 0.4f0, 0.5f0)
        modeq_f32 = ModEq(kep_f32, Float32(μ))
        @test eltype(modeq_f32) <: AbstractFloat
    end
    
    @testset "Singularity Handling" begin
        μ = 398600.4418
        
        # Test circular equatorial orbit (potential singularity)
        kep_ce = Keplerian(7000.0, 0.0, 0.0, 0.0, 0.0, 0.5)
        modeq_ce = ModEq(kep_ce, μ)
        @test all(isfinite.(params(modeq_ce)))
        @test all(.!isnan.(params(modeq_ce)))
        
        # Test equatorial orbit (i=0)
        kep_equat = Keplerian(7000.0, 0.1, 0.0, 0.0, 0.0, 0.5)
        modeq_equat = ModEq(kep_equat, μ)
        @test all(isfinite.(params(modeq_equat)))
    end
    
    @testset "Element Definitions Match Literature" begin
        μ = 398600.4418
        
        # Create known Keplerian elements
        a = 7000.0
        e = 0.1
        i = 0.2
        Ω = 0.3
        ω = 0.4
        f = 0.5
        kep = Keplerian(a, e, i, Ω, ω, f)
        
        # Convert to ModEq
        modeq = ModEq(kep, μ)
        
        # p = a(1 - e^2)
        expected_p = a * (1 - e^2)
        @test modeq.p ≈ expected_p atol=1e-8 rtol=1e-8
        
        # Verify that converting back gives same elements
        kep_back = Keplerian(modeq, μ)
        @test kep_back.a ≈ a atol=1e-6 rtol=1e-8
        @test kep_back.e ≈ e atol=1e-8 rtol=1e-8
    end
end
