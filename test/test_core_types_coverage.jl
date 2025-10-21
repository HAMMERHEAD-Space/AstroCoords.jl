using Test
using AstroCoords
using StaticArrays

@testset "Core Types Coverage" begin
    @testset "IdentityTransformation" begin
        # Test IdentityTransformation returns input unchanged
        id_trans = IdentityTransformation()
        
        # Test with Cartesian coordinate
        cart = Cartesian(1.0, 2.0, 3.0, 0.1, 0.2, 0.3)
        result = id_trans(cart)
        @test result === cart
        
        # Test with Keplerian coordinate
        kep = Keplerian(7000.0, 0.1, 0.5, 1.0, 0.5, 0.0)
        result = id_trans(kep)
        @test result === kep
        
        # Test with any object
        obj = "test_string"
        result = id_trans(obj)
        @test result === obj
    end

    @testset "ComposedTransformation" begin
        # Test composed transformation of Cartesian -> Keplerian -> Cartesian
        μ = 398600.4418
        
        # Create composition
        forward = CartesianToKeplerian
        backward = KeplerianToCartesian
        composed = ComposedTransformation(backward, forward)
        
        # Test round-trip
        cart = Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)
        result = composed(cart, μ)
        @test result isa Cartesian
        @test result.x ≈ cart.x atol=1e-8
        @test result.y ≈ cart.y atol=1e-8
        @test result.z ≈ cart.z atol=1e-8
        @test result.ẋ ≈ cart.ẋ atol=1e-8
        @test result.ẏ ≈ cart.ẏ atol=1e-8
        @test result.ż ≈ cart.ż atol=1e-8
        
        # Test composition with ∘ operator
        composed2 = backward ∘ forward
        result2 = composed2(cart, μ)
        @test result2 isa Cartesian
        @test result2.x ≈ cart.x atol=1e-8
        
        # Test three-way composition: Cart -> Kep -> ModEq -> Kep
        comp3 = KeplerianToModifiedEquinoctial ∘ CartesianToKeplerian
        result3 = comp3(cart, μ)
        @test result3 isa ModEq
    end

    @testset "Transformation Composition Rules" begin
        # Test compose function
        trans1 = CartesianToKeplerian
        trans2 = KeplerianToCartesian
        
        composed = AstroCoords.compose(trans1, trans2)
        @test composed isa ComposedTransformation
        
        # Test IdentityTransformation composition rules
        id = IdentityTransformation()
        
        # Identity composed with Identity -> Identity
        result = AstroCoords.compose(id, id)
        @test result isa IdentityTransformation
        
        # Identity composed with transformation -> transformation
        result = AstroCoords.compose(id, trans1)
        @test result === trans1
        
        result = AstroCoords.compose(trans1, id)
        @test result === trans1
    end

    @testset "Transformation Inverse" begin
        # Test inverse of ComposedTransformation
        trans1 = CartesianToKeplerian
        trans2 = KeplerianToModifiedEquinoctial
        composed = trans2 ∘ trans1
        
        inv_composed = inv(composed)
        @test inv_composed isa ComposedTransformation
        
        # Verify inverse composition order is reversed
        # inv(trans2 ∘ trans1) = inv(trans1) ∘ inv(trans2)
        # which is KeplerianToCartesian ∘ ModifiedEquinoctialToKeplerian
        
        # Test round-trip with inverse
        μ = 398600.4418
        cart = Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)
        forward = composed(cart, μ)
        backward = inv_composed(forward, μ)
        @test backward isa Cartesian
        @test backward.x ≈ cart.x atol=1e-6
        
        # Test IdentityTransformation inverse
        id = IdentityTransformation()
        inv_id = inv(id)
        @test inv_id isa IdentityTransformation
    end

    @testset "Base.zeros for Coordinate Types" begin
        # Test zeros for various Coordinate subtypes
        
        # Test with Cartesian
        z_cart = zeros(Cartesian{6, Float64})
        @test z_cart isa SVector{6, Float64}
        @test all(z_cart .== 0.0)
        
        # Test with Keplerian
        z_kep = zeros(Keplerian{6, Float64})
        @test z_kep isa SVector{6, Float64}
        @test all(z_kep .== 0.0)
        
        # Test with EDromo
        z_edromo = zeros(EDromo{8, Float64})
        @test z_edromo isa SVector{8, Float64}
        @test all(z_edromo .== 0.0)
        
        # Test with different numeric type
        z_cart32 = zeros(Cartesian{6, Float32})
        @test z_cart32 isa SVector{6, Float32}
        @test all(z_cart32 .== 0.0f0)
        
        # Test zeros with dimensions
        z_array = zeros(Cartesian{6, Float64}, 3)
        @test z_array isa Vector{SVector{6, Float64}}
        @test length(z_array) == 3
        @test all(all(z .== 0.0) for z in z_array)
        
        # Test zeros with tuple dimensions
        z_matrix = zeros(Cartesian{6, Float64}, (2, 3))
        @test z_matrix isa Matrix{SVector{6, Float64}}
        @test size(z_matrix) == (2, 3)
        
        # Test zeros with empty tuple (scalar)
        z_scalar = zeros(Cartesian{6, Float64}, ())
        @test z_scalar isa SVector{6, Float64}
    end

    @testset "AbstractTimeType Hierarchy" begin
        # Test that all time types are subtypes of AbstractTimeType
        @test PhysicalTime <: AbstractTimeType
        @test ConstantTime <: AbstractTimeType
        @test LinearTime <: AbstractTimeType
        
        # Test instantiation
        pt = PhysicalTime()
        @test pt isa PhysicalTime
        @test pt isa AbstractTimeType
        
        ct = ConstantTime()
        @test ct isa ConstantTime
        @test ct isa AbstractTimeType
        
        lt = LinearTime()
        @test lt isa LinearTime
        @test lt isa AbstractTimeType
    end

    @testset "Coordinate Type Hierarchy" begin
        # Test that coordinate types follow the correct hierarchy
        @test Cartesian <: AstroCoord
        @test Keplerian <: AstroCoord
        @test ModEq <: AstroCoord
        @test EDromo <: AstroCoord
        @test KustaanheimoStiefel <: AstroCoord
        @test StiefelScheifele <: AstroCoord
        
        # Test that AstroCoord is a Coordinate
        @test AstroCoord <: Coordinate
        
        # Test Size for Coordinate types
        @test StaticArrays.Size(Cartesian{6, Float64}) == StaticArrays.Size(6)
        @test StaticArrays.Size(EDromo{8, Float64}) == StaticArrays.Size(8)
    end

    @testset "Zero and Convert for Coordinates" begin
        # Test zero function
        cart = Cartesian(1.0, 2.0, 3.0, 0.1, 0.2, 0.3)
        z = zero(cart)
        @test z isa SVector{6, Float64}
        @test all(z .== 0.0)
        
        # Test zero with type
        z_type = zero(Cartesian{6, Float64})
        @test z_type isa SVector{6, Float64}
        @test all(z_type .== 0.0)
        
        # Test convert identity
        cart2 = Cartesian(4.0, 5.0, 6.0, 0.4, 0.5, 0.6)
        converted = convert(Cartesian{6, Float64}, cart2)
        @test converted === cart2
        
        # Test convert with different coordinate type
        kep = Keplerian(7000.0, 0.1, 0.5, 1.0, 0.5, 0.0)
        cart_from_kep = convert(Cartesian{6, Float64}, kep)
        @test cart_from_kep isa Cartesian{6, Float64}
    end
end
