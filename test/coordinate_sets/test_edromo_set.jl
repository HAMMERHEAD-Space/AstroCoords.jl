using Test
using AstroCoords
using StaticArrays

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

    @testset "Various ζ₁-ζ₈ element values" begin
        # Test with different value ranges
        
        # Positive values
        e1 = EDromo(1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.4, 100.0)
        @test all([e1.ζ₁, e1.ζ₂, e1.ζ₃, e1.ζ₄, e1.ζ₅, e1.ζ₆, e1.ζ₇, e1.ζ₈] .> 0)
        
        # Mixed positive and negative
        e2 = EDromo(-1.0, 2.0, -3.0, 0.1, -0.2, 0.3, -0.4, -100.0)
        @test e2.ζ₁ < 0
        @test e2.ζ₂ > 0
        @test e2.ζ₃ < 0
        @test e2.ζ₈ < 0
        
        # Near-zero values
        e3 = EDromo(1e-10, 2e-10, 3e-10, 1e-11, 2e-11, 3e-11, 4e-11, 1e-9)
        @test abs(e3.ζ₁) < 1e-9
        @test abs(e3.ζ₄) < 1e-10
        
        # Large values
        e4 = EDromo(1e6, 2e6, 3e6, 1e3, 2e3, 3e3, 4e3, 1e7)
        @test e4.ζ₁ > 1e5
        @test e4.ζ₈ > 1e6
    end

    @testset "Struct consistency with transformation functions" begin
        # Create a Cartesian orbit
        μ = 398600.4418  # Earth's gravitational parameter
        a = 7000.0
        e_val = 0.1
        i = 0.5
        Ω = 1.0
        ω = 0.5
        f = 0.3
        
        kep = Keplerian(a, e_val, i, Ω, ω, f)
        cart = Cartesian(kep, μ)
        
        # Create config for EDromo transformation
        config = AstroCoords.RegularizedCoordinateConfig(AstroCoords.params(cart), μ)
        ϕ = AstroCoords.compute_initial_phi(AstroCoords.params(cart), μ, config)
        
        # Convert to EDromo
        edromo_vec = AstroCoords.cart2EDromo(AstroCoords.params(cart), μ, ϕ, config)
        edromo = EDromo(edromo_vec...)
        
        # Verify struct fields match vector elements
        @test edromo.ζ₁ ≈ edromo_vec[1]
        @test edromo.ζ₂ ≈ edromo_vec[2]
        @test edromo.ζ₃ ≈ edromo_vec[3]
        @test edromo.ζ₄ ≈ edromo_vec[4]
        @test edromo.ζ₅ ≈ edromo_vec[5]
        @test edromo.ζ₆ ≈ edromo_vec[6]
        @test edromo.ζ₇ ≈ edromo_vec[7]
        @test edromo.ζ₈ ≈ edromo_vec[8]
        
        # Verify params() extracts all fields correctly
        params_vec = AstroCoords.params(edromo)
        @test params_vec[1] ≈ edromo.ζ₁
        @test params_vec[2] ≈ edromo.ζ₂
        @test params_vec[3] ≈ edromo.ζ₃
        @test params_vec[4] ≈ edromo.ζ₄
        @test params_vec[5] ≈ edromo.ζ₅
        @test params_vec[6] ≈ edromo.ζ₆
        @test params_vec[7] ≈ edromo.ζ₇
        @test params_vec[8] ≈ edromo.ζ₈
        
        # Verify round-trip: EDromo → Cartesian → EDromo
        cart_recovered_vec = AstroCoords.EDromo2cart(AstroCoords.params(edromo), μ, ϕ, config)
        cart_recovered = Cartesian(cart_recovered_vec...)
        
        @test cart_recovered.x ≈ cart.x rtol=1e-9
        @test cart_recovered.y ≈ cart.y rtol=1e-9
        @test cart_recovered.z ≈ cart.z rtol=1e-9
        @test cart_recovered.ẋ ≈ cart.ẋ rtol=1e-9
        @test cart_recovered.ẏ ≈ cart.ẏ rtol=1e-9
        @test cart_recovered.ż ≈ cart.ż rtol=1e-9
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
