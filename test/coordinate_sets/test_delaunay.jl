using Test
using AstroCoords
using StaticArrays

@testset "Delaunay Coordinate Set" begin
    @testset "Inner constructor Delaunay{T}" begin
        d = Delaunay{Float64}(1.0, 2.0, 3.0, 0.5, 1.0, 2.0)
        @test d isa Delaunay{Float64}
        @test d.L ≈ 1.0
        @test d.G ≈ 2.0
        @test d.H ≈ 3.0
        @test d.M ≈ 0.5
        @test d.ω ≈ 1.0
        @test d.Ω ≈ 2.0
    end

    @testset "Constructor with type promotion" begin
        # Float32 + Float64 should promote to Float64
        d = Delaunay(1.0f0, 2.0, 3.0f0, 0.5, 1.0f0, 2.0)
        @test d isa Delaunay{Float64}
        @test d.L ≈ 1.0
        @test d.G ≈ 2.0
        @test d.H ≈ 3.0
        @test d.M ≈ 0.5
        @test d.ω ≈ 1.0
        @test d.Ω ≈ 2.0
    end

    @testset "StaticVector constructor" begin
        vec = SVector{6}(1.0, 2.0, 3.0, 0.5, 1.0, 2.0)
        d = Delaunay(vec)
        @test d isa Delaunay{Float64}
        @test d.L ≈ 1.0
        @test d.G ≈ 2.0
        @test d.H ≈ 3.0
        @test d.M ≈ 0.5
        @test d.ω ≈ 1.0
        @test d.Ω ≈ 2.0
    end

    @testset "Constructor from vector" begin
        vec = [1.0, 2.0, 3.0, 0.5, 1.0, 2.0]
        d = Delaunay(vec)
        @test d isa Delaunay{Float64}
        @test d.L ≈ 1.0
        @test d.G ≈ 2.0
        @test d.H ≈ 3.0
        @test d.M ≈ 0.5
        @test d.ω ≈ 1.0
        @test d.Ω ≈ 2.0
    end

    @testset "Base.one(Delaunay)" begin
        d_one = Base.one(Delaunay)
        @test d_one isa Delaunay{Float64}
        @test d_one.L ≈ 0.0
        @test d_one.G ≈ 0.0
        @test d_one.H ≈ 0.0
        @test d_one.M ≈ 0.0
        @test d_one.ω ≈ 0.0
        @test d_one.Ω ≈ 0.0

        d_one_f32 = Base.one(Delaunay; T=Float32)
        @test d_one_f32 isa Delaunay{Float32}
        @test d_one_f32.L ≈ 0.0f0
    end

    @testset "Base.getindex for all 6 indices" begin
        d = Delaunay(1.0, 2.0, 3.0, 0.5, 1.0, 2.0)
        
        @test d[1] ≈ 1.0  # L
        @test d[2] ≈ 2.0  # G
        @test d[3] ≈ 3.0  # H
        @test d[4] ≈ 0.5  # M
        @test d[5] ≈ 1.0  # ω
        @test d[6] ≈ 2.0  # Ω
        
        @test_throws BoundsError d[0]
        @test_throws BoundsError d[7]
    end

    @testset "Canonical Hamiltonian properties" begin
        # Delaunay elements are canonical (action-angle) variables
        # L = √(μa), G = L√(1-e²), H = G cos(i)
        μ = 398600.4418  # Earth's gravitational parameter (km³/s²)
        
        # Example elliptical orbit
        a = 7000.0  # km
        e = 0.1
        i = 0.5  # rad
        Ω = 1.0  # rad
        ω = 0.5  # rad
        M = 0.3  # rad
        
        L_expected = √(μ * a)
        G_expected = L_expected * √(1 - e^2)
        H_expected = G_expected * cos(i)
        
        d = Delaunay(L_expected, G_expected, H_expected, M, ω, Ω)
        
        @test d.L ≈ L_expected rtol=1e-10
        @test d.G ≈ G_expected rtol=1e-10
        @test d.H ≈ H_expected rtol=1e-10
        @test d.M ≈ M rtol=1e-10
        @test d.ω ≈ ω rtol=1e-10
        @test d.Ω ≈ Ω rtol=1e-10
    end

    @testset "Various orbit types" begin
        μ = 398600.4418
        
        # Circular orbit (e=0)
        a = 7000.0
        e = 0.0
        i = 0.5
        L = √(μ * a)
        G = L * √(1 - e^2)
        H = G * cos(i)
        d_circ = Delaunay(L, G, H, 0.0, 0.0, 0.0)
        @test d_circ.L ≈ d_circ.G rtol=1e-10  # For circular orbit, L=G
        
        # Equatorial orbit (i=0)
        i_eq = 0.0
        H_eq = G * cos(i_eq)
        d_eq = Delaunay(L, G, H_eq, 0.0, 0.0, 0.0)
        @test d_eq.H ≈ d_eq.G rtol=1e-10  # For equatorial orbit, H=G
        
        # Highly eccentric orbit
        e_high = 0.9
        G_high = L * √(1 - e_high^2)
        H_high = G_high * cos(i)
        d_ecc = Delaunay(L, G_high, H_high, 0.0, 0.0, 0.0)
        @test d_ecc.G < d_ecc.L  # High eccentricity means G << L
    end

    @testset "Action-angle variable definitions" begin
        μ = 398600.4418
        
        # Verify action variables (L, G, H) are related by orbital parameters
        a = 8000.0
        e = 0.2
        i = 0.7
        
        L = √(μ * a)
        G = L * √(1 - e^2)
        H = G * cos(i)
        
        # Verify inequalities that must hold
        @test H ≤ G  # H ≤ G always (since |cos(i)| ≤ 1)
        @test G ≤ L  # G ≤ L always (since √(1-e²) ≤ 1)
        @test H ≤ L  # Transitive property
        
        # Verify eccentricity recovery
        e_recovered = √(1 - (G/L)^2)
        @test e_recovered ≈ e rtol=1e-10
        
        # Verify inclination recovery
        i_recovered = acos(H/G)
        @test i_recovered ≈ i rtol=1e-10
    end
end
