using Test
using AstroCoords
using StaticArrays

@testset "Core Type System" begin
    @testset "ComposedTransformation Construction" begin
        # Test basic composition (line 54)
        trans1 = AstroCoords.CartesianToKeplerian
        trans2 = AstroCoords.KeplerianToModEq
        composed = AstroCoords.ComposedTransformation(trans2, trans1)

        @test composed isa AstroCoords.ComposedTransformation
        @test composed.t1 === trans2
        @test composed.t2 === trans1
    end

    @testset "ComposedTransformation Application" begin
        # Test transformation application (lines 58-59)
        μ = 398600.4418  # Earth μ (km³/s²)
        cart = Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)

        # Compose Cart→Kep→ModEq
        trans = AstroCoords.ComposedTransformation(
            AstroCoords.KeplerianToCartesian, AstroCoords.CartesianToKeplerian
        )
        result = trans(cart, μ)

        @test result isa Cartesian
        # Should be round-trip (within tolerance)
        @test result.x ≈ cart.x atol=1e-6
        @test result.y ≈ cart.y atol=1e-6
        @test result.z ≈ cart.z atol=1e-6
    end

    @testset "Transformation Composition Chains" begin
        # Test composition operator (lines 66-67)
        trans1 = AstroCoords.CartesianToKeplerian
        trans2 = AstroCoords.KeplerianToModEq
        trans3 = AstroCoords.ModEqToKeplerian

        # Test ∘ operator
        composed_op = trans1 ∘ trans2
        @test composed_op isa AstroCoords.ComposedTransformation

        # Test compose function
        composed_fn = AstroCoords.compose(trans1, trans2)
        @test composed_fn isa AstroCoords.ComposedTransformation

        # Test chaining multiple compositions
        chain = trans1 ∘ trans2 ∘ trans3
        @test chain isa AstroCoords.ComposedTransformation
        @test chain.t1 isa
            Union{AstroCoords.Transformation,AstroCoords.ComposedTransformation}
    end

    @testset "Transformation with Various Coordinate Types" begin
        # Test type combinations (lines 82-90)
        μ = 398600.4418

        # Test with different coordinate types
        cart = Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)
        kep = Keplerian(cart, μ)
        modeq = ModEq(kep, μ)

        # Compose transformations between different types
        cart_to_modeq = AstroCoords.KeplerianToModEq ∘ AstroCoords.CartesianToKeplerian
        result = cart_to_modeq(cart, μ)
        @test result isa ModEq

        # Test round-trip with composition
        modeq_to_cart = AstroCoords.KeplerianToCartesian ∘ AstroCoords.ModEqToKeplerian
        back = modeq_to_cart(result, μ)
        @test back isa Cartesian
        @test back.x ≈ cart.x atol=1e-6
    end

    @testset "IdentityTransformation" begin
        # Test IdentityTransformation behavior
        id = AstroCoords.IdentityTransformation()

        cart = Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)
        result = id(cart)

        @test result === cart
        @test id(42) === 42
        @test id("test") === "test"
    end

    @testset "Transformation Type Stability" begin
        # Test type stability of composed transformations
        μ = 398600.4418
        cart = Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)

        trans = AstroCoords.KeplerianToCartesian ∘ AstroCoords.CartesianToKeplerian
        result = trans(cart, μ)

        @test typeof(result) == Cartesian{Float64}
        @test eltype(params(result)) == Float64
    end

    @testset "Transformation Algebra - Associativity" begin
        # Test composition associativity: (A ∘ B) ∘ C == A ∘ (B ∘ C)
        # For transformations: Cart → Kep → ModEq → Kep
        # We compose as: t3 ∘ t2 ∘ t1 (right to left application)
        μ = 398600.4418
        cart = Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)

        t1 = AstroCoords.CartesianToKeplerian
        t2 = AstroCoords.KeplerianToModEq
        t3 = AstroCoords.ModEqToKeplerian

        # Left association: (t3 ∘ t2) ∘ t1
        left = (t3 ∘ t2) ∘ t1
        result_left = left(cart, μ)

        # Right association: t3 ∘ (t2 ∘ t1)
        right = t3 ∘ (t2 ∘ t1)
        result_right = right(cart, μ)

        @test result_left isa Keplerian
        @test result_right isa Keplerian
        @test result_left.a ≈ result_right.a atol=1e-10
        @test result_left.e ≈ result_right.e atol=1e-10
    end

    @testset "Transformation Inverse" begin
        # Test inverse transformations
        trans = AstroCoords.CartesianToKeplerian
        inv_trans = inv(trans)

        @test inv_trans === AstroCoords.KeplerianToCartesian
        @test inv(inv_trans) === trans

        # Test composed transformation inverse
        composed = AstroCoords.CartesianToKeplerian ∘ AstroCoords.KeplerianToModEq
        inv_composed = inv(composed)
        @test inv_composed isa AstroCoords.ComposedTransformation
    end

    @testset "IdentityTransformation Composition" begin
        # Test identity composition rules
        id = AstroCoords.IdentityTransformation()
        trans = AstroCoords.CartesianToKeplerian

        # compose(id, id) = id
        @test AstroCoords.compose(id, id) === id

        # compose(id, trans) = trans
        @test AstroCoords.compose(id, trans) === trans

        # compose(trans, id) = trans
        @test AstroCoords.compose(trans, id) === trans
    end
end
