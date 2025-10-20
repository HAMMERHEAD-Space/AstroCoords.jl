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
        k = Keplerian(7000.0f0, 0.01, 0.1f0, 0.2, 0.5f0, 0.4)
        @test k isa Keplerian{Float64}
        @test k.a ≈ 7000.0
        @test k.e ≈ 0.01
        @test k.i ≈ 0.1
        @test k.Ω ≈ 0.2
        @test k.ω ≈ 0.5
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

    @testset "Type stability" begin
        k_f64 = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.4)
        @test typeof(k_f64.a) === Float64
        @test typeof(k_f64.e) === Float64

        k_f32 = Keplerian(7000.0f0, 0.01f0, 0.1f0, 0.2f0, 0.3f0, 0.4f0)
        @test typeof(k_f32.a) === Float32
        @test typeof(k_f32.e) === Float32
    end

    @testset "Exotic coordinate properties" begin
        @testset "Basic property access" begin
            k = Keplerian(7000.0, 0.1, 0.0, 0.0, 0.0, π/4)

            # Test that regular fields still work
            @test k.a ≈ 7000.0
            @test k.e ≈ 0.1
            @test k.f ≈ π/4

            # Test that exotic properties are not fields but are accessible
            @test hasfield(typeof(k), :M) == false
            @test hasfield(typeof(k), :E) == false
            @test hasmethod(getproperty, (typeof(k), Symbol))

            # Test property values are computed correctly
            M_expected = trueAnomaly2MeanAnomaly(k.f, k.e)
            E_expected = trueAnomaly2EccentricAnomaly(k.f, k.e)
            @test k.M ≈ M_expected
            @test k.E ≈ E_expected
        end

        @testset "Property names" begin
            k = Keplerian(7000.0, 0.1, 0.0, 0.0, 0.0, π/4)
            props = propertynames(k)

            # Test that all expected properties are included
            @test :a in props
            @test :e in props
            @test :i in props
            @test :Ω in props
            @test :ω in props
            @test :f in props
            @test :M in props  # Mean anomaly
            @test :E in props  # Eccentric anomaly

            # Test exact property list
            @test props == (:a, :e, :i, :Ω, :ω, :f, :M, :E)
        end

        @testset "Round-trip consistency" begin
            k = Keplerian(7000.0, 0.2, 0.1, 0.2, 0.3, π/3)

            # Test that we can round-trip through mean anomaly
            M = k.M
            f_roundtrip = meanAnomaly2TrueAnomaly(M, k.e)
            @test f_roundtrip ≈ k.f atol=1e-12

            # Test that we can round-trip through eccentric anomaly
            E = k.E
            f_roundtrip2 = eccentricAnomaly2TrueAnomaly(E, k.e)
            @test f_roundtrip2 ≈ k.f atol=1e-12
        end
    end
end
