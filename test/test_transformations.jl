@testset "Transformation Infrastructure" begin
    μ = 398600.4418  # Earth standard gravitational parameter

    @testset "Transformation Type Hierarchy" begin
        # Verify transformations are subtypes of correct abstract types
        @test AstroCoords.CartesianToKeplerianTransform <:
            AstroCoords.AstroCoordTransformation
        @test AstroCoords.AstroCoordTransformation <:
            AstroCoords.AstrodynamicsTransformation
        @test AstroCoords.AstrodynamicsTransformation <: AstroCoords.Transformation
    end

    @testset "Transformation Inverse Relationship" begin
        # Test that inv() returns the correct inverse transformation
        cart_to_kep = AstroCoords.CartesianToKeplerian
        kep_to_cart = AstroCoords.KeplerianToCartesian

        @test inv(cart_to_kep) === kep_to_cart
        @test inv(kep_to_cart) === cart_to_kep
    end

    @testset "IdentityTransformation" begin
        # Test identity transformation
        id_trans = AstroCoords.IdentityTransformation()

        # Create test coordinate
        cart = Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)

        # Identity should return input unchanged
        result = id_trans(cart)
        @test result === cart

        # Inverse of identity is identity
        @test inv(id_trans) === id_trans
    end

    @testset "ComposedTransformation - Cartesian → Keplerian → ModEq" begin
        # Test composed transformation
        cart_to_kep = AstroCoords.CartesianToKeplerian
        kep_to_modeq = AstroCoords.KeplerianToModEq

        composed = kep_to_modeq ∘ cart_to_kep
        @test composed isa AstroCoords.ComposedTransformation

        # Test on actual data
        cart = Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)
        result = composed(cart, μ)

        @test result isa ModEq
    end

    @testset "ComposedTransformation - Round Trip" begin
        # Test Cartesian → Keplerian → Cartesian
        cart_to_kep = AstroCoords.CartesianToKeplerian
        kep_to_cart = AstroCoords.KeplerianToCartesian

        round_trip = kep_to_cart ∘ cart_to_kep

        cart_orig = Cartesian(7000.0, 1000.0, 500.0, 0.5, 7.0, 0.3)
        cart_back = round_trip(cart_orig, μ)

        # Should recover original (within tolerance)
        @test cart_back.x ≈ cart_orig.x atol=1e-9
        @test cart_back.y ≈ cart_orig.y atol=1e-9
        @test cart_back.z ≈ cart_orig.z atol=1e-9
        @test cart_back.ẋ ≈ cart_orig.ẋ atol=1e-9
        @test cart_back.ẏ ≈ cart_orig.ẏ atol=1e-9
        @test cart_back.ż ≈ cart_orig.ż atol=1e-9
    end

    @testset "ComposedTransformation - Inverse" begin
        # Test that inverse of composed is reverse composition
        trans1 = AstroCoords.CartesianToKeplerian
        trans2 = AstroCoords.KeplerianToModEq

        composed = trans1 ∘ trans2
        inv_composed = inv(composed)

        # inv((A ∘ B)) = inv(B) ∘ inv(A)
        @test inv_composed isa AstroCoords.ComposedTransformation
    end

    @testset "Transformation Type Stability" begin
        # Test that transformations maintain type stability
        cart = Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)

        # Cartesian to Keplerian
        kep = AstroCoords.CartesianToKeplerian(cart, μ)
        @test kep isa Keplerian{Float64}

        # Back to Cartesian
        cart_back = AstroCoords.KeplerianToCartesian(kep, μ)
        @test cart_back isa Cartesian{Float64}
    end

    @testset "Transformation Type Promotion (Float32/Float64)" begin
        # Test type promotion with mixed Float32/Float64
        cart = Cartesian(
            Float32(7000.0),
            Float32(0.0),
            Float32(0.0),
            Float32(0.0),
            Float32(7.5),
            Float32(0.0),
        )
        μ_f64 = Float64(398600.4418)

        kep = AstroCoords.CartesianToKeplerian(cart, μ_f64)
        @test kep isa Keplerian{Float64}  # Should promote to Float64
    end

    @testset "Direct Transformation Pairs Exist" begin
        # Verify all defined transformation pairs exist
        @test isdefined(AstroCoords, :CartesianToKeplerian)
        @test isdefined(AstroCoords, :KeplerianToCartesian)
        @test isdefined(AstroCoords, :CartesianToMilankovich)
        @test isdefined(AstroCoords, :MilankovichToCartesian)
        @test isdefined(AstroCoords, :CartesianToCylindrical)
        @test isdefined(AstroCoords, :CylindricalToCartesian)
        @test isdefined(AstroCoords, :CartesianToSpherical)
        @test isdefined(AstroCoords, :SphericalToCartesian)
        @test isdefined(AstroCoords, :CartesianToDelaunay)
        @test isdefined(AstroCoords, :DelaunayToCartesian)
        @test isdefined(AstroCoords, :KeplerianToUSM7)
        @test isdefined(AstroCoords, :USM7ToKeplerian)
        @test isdefined(AstroCoords, :KeplerianToModEq)
        @test isdefined(AstroCoords, :ModEqToKeplerian)
    end

    @testset "Transformation Aliases" begin
        # Test that readable aliases exist
        @test isdefined(AstroCoords, :KeplerianToModifiedEquinoctial)
        @test isdefined(AstroCoords, :ModifiedEquinoctialToKeplerian)

        # Verify aliases point to same transformations
        @test AstroCoords.KeplerianToModifiedEquinoctial === AstroCoords.KeplerianToModEq
        @test AstroCoords.ModifiedEquinoctialToKeplerian === AstroCoords.ModEqToKeplerian
    end

    @testset "Regularized Coordinate Transformations" begin
        # Test that regularized coordinate transformations exist
        @test isdefined(AstroCoords, :CartesianToEDromoTransform)
        @test isdefined(AstroCoords, :EDromoToCartesianTransform)
        @test isdefined(AstroCoords, :CartesianToKustaanheimoStiefelTransform)
        @test isdefined(AstroCoords, :KustaanheimoStiefelToCartesianTransform)
        @test isdefined(AstroCoords, :CartesianToStiefelScheifeleTransform)
        @test isdefined(AstroCoords, :StiefelScheifeleToCartesianTransform)
    end

    @testset "Transformation Call Signature" begin
        # Test that transformations accept correct arguments
        cart = Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)
        trans = AstroCoords.CartesianToKeplerian

        # Should accept (coordinate, μ)
        result = trans(cart, μ)
        @test result isa Keplerian

        # Test with additional args (for regularized coords)
        config = RegularizedCoordinateConfig(DU=1.0, TU=1.0, W=0.0, t₀=0.0)
        # This should work for compatible transformations
        # (not all transformations need config, but they should accept it)
    end

    @testset "Transformation Error Handling" begin
        # Test that invalid transformations produce sensible errors
        trans = AstroCoords.IdentityTransformation()

        # inv() should be defined for all transformations
        @test inv(trans) isa AstroCoords.Transformation
    end

    @testset "Multiple Composition" begin
        # Test composing more than two transformations
        trans1 = AstroCoords.CartesianToKeplerian
        trans2 = AstroCoords.KeplerianToModEq
        trans3 = AstroCoords.ModEqToKeplerian

        composed = trans3 ∘ trans2 ∘ trans1
        @test composed isa AstroCoords.ComposedTransformation

        # Test on data
        cart = Cartesian(7000.0, 1000.0, 500.0, 0.5, 7.0, 0.3)
        result = composed(cart, μ)
        @test result isa Keplerian
    end

    @testset "Identity Composition Rules" begin
        # Test that identity composes correctly
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
