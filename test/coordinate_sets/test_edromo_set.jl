@testset "EDromo Coordinate Set" begin
    @testset "Inner constructor EDromo{T}" begin
        e = EDromo{Float64}(1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.4, 100.0)
        @test e isa EDromo{Float64}
        @test e.ζ₁ ≈ 1.0
        @test e.ζ₂ ≈ 2.0
        @test e.ζ₃ ≈ 3.0
        @test e.ζ₄ ≈ 0.1
        @test e.ζ₅ ≈ 0.2
        @test e.ζ₆ ≈ 0.3
        @test e.ζ₇ ≈ 0.4
        @test e.ζ₈ ≈ 100.0
    end

    @testset "Constructor from AbstractVector" begin
        vec = [1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.4, 100.0]
        e = EDromo(vec)
        @test e isa EDromo{Float64}
        @test e.ζ₁ ≈ 1.0
        @test e.ζ₂ ≈ 2.0
        @test e.ζ₃ ≈ 3.0
        @test e.ζ₄ ≈ 0.1
        @test e.ζ₅ ≈ 0.2
        @test e.ζ₆ ≈ 0.3
        @test e.ζ₇ ≈ 0.4
        @test e.ζ₈ ≈ 100.0
    end

    @testset "Constructor with type promotion" begin
        # Mix Float32 and Float64 - should promote to Float64
        e = EDromo(1.0f0, 2.0, 3.0f0, 0.1, 0.2f0, 0.3, 0.4f0, 100.0)
        @test e isa EDromo{Float64}
        @test typeof(e.ζ₁) === Float64
        @test typeof(e.ζ₂) === Float64
        @test e.ζ₁ ≈ 1.0
        @test e.ζ₂ ≈ 2.0
        @test e.ζ₃ ≈ 3.0
        @test e.ζ₄ ≈ 0.1
        @test e.ζ₅ ≈ 0.2
        @test e.ζ₆ ≈ 0.3
        @test e.ζ₇ ≈ 0.4
        @test e.ζ₈ ≈ 100.0
    end

    @testset "Constructor with individual parameters" begin
        e = EDromo(1.5, 2.5, 3.5, 0.15, 0.25, 0.35, 0.45, 150.0)
        @test e isa EDromo{Float64}
        @test e.ζ₁ ≈ 1.5
        @test e.ζ₂ ≈ 2.5
        @test e.ζ₃ ≈ 3.5
        @test e.ζ₄ ≈ 0.15
        @test e.ζ₅ ≈ 0.25
        @test e.ζ₆ ≈ 0.35
        @test e.ζ₇ ≈ 0.45
        @test e.ζ₈ ≈ 150.0
    end

    @testset "StaticVector constructor" begin
        vec = SVector{8}(1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.4, 100.0)
        e = EDromo(vec)
        @test e isa EDromo{Float64}
        @test e.ζ₁ ≈ 1.0
        @test e.ζ₂ ≈ 2.0
        @test e.ζ₃ ≈ 3.0
        @test e.ζ₄ ≈ 0.1
        @test e.ζ₅ ≈ 0.2
        @test e.ζ₆ ≈ 0.3
        @test e.ζ₇ ≈ 0.4
        @test e.ζ₈ ≈ 100.0
    end

    @testset "Base.getindex for all 8 indices" begin
        e = EDromo(1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.4, 100.0)

        @test e[1] ≈ 1.0   # ζ₁
        @test e[2] ≈ 2.0   # ζ₂
        @test e[3] ≈ 3.0   # ζ₃
        @test e[4] ≈ 0.1   # ζ₄
        @test e[5] ≈ 0.2   # ζ₅
        @test e[6] ≈ 0.3   # ζ₆
        @test e[7] ≈ 0.4   # ζ₇
        @test e[8] ≈ 100.0 # ζ₈

        @test_throws BoundsError e[0]
        @test_throws BoundsError e[9]
    end

    @testset "Base.one(EDromo)" begin
        e_one = Base.one(EDromo)
        @test e_one isa EDromo{Float64}
        @test e_one.ζ₁ ≈ 0.0
        @test e_one.ζ₂ ≈ 0.0
        @test e_one.ζ₃ ≈ 0.0
        @test e_one.ζ₄ ≈ 0.0
        @test e_one.ζ₅ ≈ 0.0
        @test e_one.ζ₆ ≈ 0.0
        @test e_one.ζ₇ ≈ 0.0
        @test e_one.ζ₈ ≈ 0.0

        e_one_f32 = Base.one(EDromo; T=Float32)
        @test e_one_f32 isa EDromo{Float32}
        @test e_one_f32.ζ₁ ≈ 0.0f0
    end

    @testset "Type stability" begin
        e_f64 = EDromo(1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.4, 100.0)
        @test typeof(e_f64.ζ₁) === Float64
        @test typeof(e_f64.ζ₈) === Float64

        e_f32 = EDromo(1.0f0, 2.0f0, 3.0f0, 0.1f0, 0.2f0, 0.3f0, 0.4f0, 100.0f0)
        @test typeof(e_f32.ζ₁) === Float32
        @test typeof(e_f32.ζ₈) === Float32
    end
end
