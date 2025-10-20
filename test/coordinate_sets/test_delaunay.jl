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

    @testset "Exotic coordinate properties" begin
        @testset "Basic property access" begin
            # Create realistic Delaunay coordinates
            μ = 3.986004418e14  # Earth's gravitational parameter m³/s²
            a = 7000e3          # Semi-major axis in meters
            e = 0.2             # Eccentricity
            i = 0.1             # Inclination
            L = √(μ * a)
            G = √(μ * a * (1 - e^2))
            H = G * cos(i)
            M = π/4             # Mean anomaly
            ω = 0.2             # Argument of periapsis
            Ω = 0.3             # RAAN

            del = Delaunay(L, G, H, M, ω, Ω)

            # Test that regular fields still work
            @test del.L ≈ L
            @test del.G ≈ G
            @test del.M ≈ M

            # Test that exotic properties are not fields but are accessible
            @test hasfield(typeof(del), :E) == false
            @test hasfield(typeof(del), :f) == false
            @test hasmethod(getproperty, (typeof(del), Symbol))

            # Test property values are computed correctly
            e_computed = √(1 - (G/L)^2)
            E_expected = meanAnomaly2EccentricAnomaly(M, e_computed)
            f_expected = meanAnomaly2TrueAnomaly(M, e_computed)
            @test del.E ≈ E_expected
            @test del.f ≈ f_expected

            # Verify eccentricity extraction is correct
            @test e_computed ≈ e atol=1e-14
        end

        @testset "Property names" begin
            del = Delaunay(1.0, 0.9, 0.8, π/6, 0.0, 0.0)
            props = propertynames(del)

            # Test that all expected properties are included
            @test :L in props
            @test :G in props
            @test :H in props
            @test :M in props
            @test :ω in props
            @test :Ω in props
            @test :E in props  # Eccentric anomaly
            @test :f in props  # True anomaly

            # Test exact property list
            @test props == (:L, :G, :H, :M, :ω, :Ω, :E, :f)
        end

        @testset "Round-trip consistency" begin
            # Create Delaunay coordinates
            L = 2.0
            G = 1.8
            H = 1.6
            M = π/3
            del = Delaunay(L, G, H, M, 0.1, 0.2)

            # Test that we can round-trip through the anomalies
            e_computed = √(1 - (G/L)^2)
            E = del.E
            f = del.f

            # Round-trip E -> f -> E
            f_from_E = eccentricAnomaly2TrueAnomaly(E, e_computed)
            E_from_f = trueAnomaly2EccentricAnomaly(f_from_E, e_computed)
            @test E ≈ E_from_f atol=1e-12

            # Round-trip M -> f -> M
            f_from_M = meanAnomaly2TrueAnomaly(M, e_computed)
            M_from_f = trueAnomaly2MeanAnomaly(f_from_M, e_computed)
            @test M ≈ M_from_f atol=1e-12
        end
    end
end
