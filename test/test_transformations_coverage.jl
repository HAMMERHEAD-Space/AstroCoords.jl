using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "@define_transformation_pair macro instantiation" begin
    # Test that the macro generates working transformations
    # We'll test Cartesian ↔ Keplerian as an example
    
    μ = 398600.4418
    r = [7000.0, 1000.0, 2000.0]
    v = [1.0, 7.0, -0.5]
    cart = Cartesian(vcat(r, v)...)
    
    # Forward transformation
    kep = CartesianToKeplerian(cart, μ)
    @test kep isa Keplerian
    
    # Inverse transformation
    cart_back = KeplerianToCartesian(kep, μ)
    @test cart_back isa Cartesian
    
    # Round-trip should preserve values
    @test params(cart_back) ≈ params(cart) atol=1e-8 rtol=1e-10
end

@testset "Transformation objects: Cartesian ↔ Keplerian" begin
    μ = 398600.4418
    cart_coord = Cartesian(7000.0, 1000.0, 2000.0, 1.0, 7.0, -0.5)
    
    # Test the transformation object
    transform = CartesianToKeplerian
    @test transform isa AstroCoordTransformation
    
    # Apply transformation
    kep_coord = transform(cart_coord, μ)
    @test kep_coord isa Keplerian
    
    # Test inverse
    inv_transform = inv(transform)
    @test inv_transform isa AstroCoordTransformation
    cart_back = inv_transform(kep_coord, μ)
    @test params(cart_back) ≈ params(cart_coord) atol=1e-8 rtol=1e-10
end

@testset "Transformation objects: Cartesian ↔ ModEq" begin
    μ = 398600.4418
    cart_coord = Cartesian(7000.0, 1000.0, 2000.0, 1.0, 7.0, -0.5)
    
    # Test Cartesian → Keplerian → ModEq path
    kep_coord = CartesianToKeplerian(cart_coord, μ)
    modeq_coord = KeplerianToModEq(kep_coord, μ)
    @test modeq_coord isa ModEq
    
    # Test inverse path
    kep_back = ModEqToKeplerian(modeq_coord, μ)
    cart_back = KeplerianToCartesian(kep_back, μ)
    @test params(cart_back) ≈ params(cart_coord) atol=1e-6 rtol=1e-8
end

@testset "Transformation objects: Cartesian ↔ Spherical" begin
    μ = 398600.4418
    cart_coord = Cartesian(7000.0, 1000.0, 2000.0, 1.0, 7.0, -0.5)
    
    sphere_coord = CartesianToSpherical(cart_coord, μ)
    @test sphere_coord isa Spherical
    
    cart_back = SphericalToCartesian(sphere_coord, μ)
    @test params(cart_back) ≈ params(cart_coord) atol=1e-8 rtol=1e-10
end

@testset "Transformation objects: Cartesian ↔ Cylindrical" begin
    μ = 398600.4418
    cart_coord = Cartesian(7000.0, 1000.0, 2000.0, 1.0, 7.0, -0.5)
    
    cylind_coord = CartesianToCylindrical(cart_coord, μ)
    @test cylind_coord isa Cylindrical
    
    cart_back = CylindricalToCartesian(cylind_coord, μ)
    @test params(cart_back) ≈ params(cart_coord) atol=1e-8 rtol=1e-10
end

@testset "Transformation objects: Cartesian ↔ Delaunay" begin
    μ = 398600.4418
    cart_coord = Cartesian(7000.0, 1000.0, 2000.0, 1.0, 7.0, -0.5)
    
    del_coord = CartesianToDelaunay(cart_coord, μ)
    @test del_coord isa Delaunay
    
    cart_back = DelaunayToCartesian(del_coord, μ)
    @test params(cart_back) ≈ params(cart_coord) atol=1e-6 rtol=1e-8
end

@testset "Transformation objects: Cartesian ↔ Milankovich" begin
    μ = 398600.4418
    cart_coord = Cartesian(7000.0, 1000.0, 2000.0, 1.0, 7.0, -0.5)
    
    mil_coord = CartesianToMilankovich(cart_coord, μ)
    @test mil_coord isa Milankovich
    
    cart_back = MilankovichToCartesian(mil_coord, μ)
    @test params(cart_back) ≈ params(cart_coord) atol=1e-6 rtol=1e-8
end

@testset "Transformation objects: Keplerian ↔ USM7" begin
    μ = 398600.4418
    kep_coord = Keplerian(7000.0, 0.2, 0.5, 1.0, 0.5, 1.0)
    
    usm7_coord = KeplerianToUSM7(kep_coord, μ)
    @test usm7_coord isa USM7
    
    kep_back = USM7ToKeplerian(usm7_coord, μ)
    @test params(kep_back) ≈ params(kep_coord) atol=1e-8 rtol=1e-10
end

@testset "Transformation objects: USM7 ↔ USM6" begin
    μ = 398600.4418
    kep_coord = Keplerian(7000.0, 0.2, 0.5, 1.0, 0.5, 1.0)
    usm7_coord = KeplerianToUSM7(kep_coord, μ)
    
    usm6_coord = USM7ToUSM6(usm7_coord, μ)
    @test usm6_coord isa USM6
    
    usm7_back = USM6ToUSM7(usm6_coord, μ)
    @test params(usm7_back) ≈ params(usm7_coord) atol=1e-10 rtol=1e-10
end

@testset "Transformation objects: USM7 ↔ USMEM" begin
    μ = 398600.4418
    kep_coord = Keplerian(7000.0, 0.2, 0.5, 1.0, 0.5, 1.0)
    usm7_coord = KeplerianToUSM7(kep_coord, μ)
    
    usmem_coord = USM7ToUSMEM(usm7_coord, μ)
    @test usmem_coord isa USMEM
    
    usm7_back = USMEMToUSM7(usmem_coord, μ)
    @test params(usm7_back) ≈ params(usm7_coord) atol=1e-10 rtol=1e-10
end

@testset "IdentityTransformation behavior" begin
    μ = 398600.4418
    cart_coord = Cartesian(7000.0, 1000.0, 2000.0, 1.0, 7.0, -0.5)
    
    # IdentityTransformation should return input unchanged
    id_trans = IdentityTransformation()
    result = id_trans(cart_coord)
    @test result === cart_coord
    
    # inv(IdentityTransformation) should be IdentityTransformation
    @test inv(id_trans) isa IdentityTransformation
end

@testset "ComposedTransformation: Cartesian → Keplerian → ModEq" begin
    μ = 398600.4418
    cart_coord = Cartesian(7000.0, 1000.0, 2000.0, 1.0, 7.0, -0.5)
    
    # Compose transformations
    composed = CartesianToKeplerian ∘ KeplerianToModEq
    @test composed isa ComposedTransformation
    
    # This should fail because the order is wrong - we need to apply in reverse
    # The correct composition should be: KeplerianToModEq ∘ CartesianToKeplerian
    composed_correct = KeplerianToModEq ∘ CartesianToKeplerian
    
    # Apply composed transformation
    modeq_coord = composed_correct(cart_coord, μ)
    @test modeq_coord isa ModEq
    
    # Verify result matches step-by-step application
    kep_coord = CartesianToKeplerian(cart_coord, μ)
    modeq_direct = KeplerianToModEq(kep_coord, μ)
    @test params(modeq_coord) ≈ params(modeq_direct) atol=1e-12 rtol=1e-12
end

@testset "ComposedTransformation: inverse" begin
    μ = 398600.4418
    cart_coord = Cartesian(7000.0, 1000.0, 2000.0, 1.0, 7.0, -0.5)
    
    # Compose transformations
    composed = KeplerianToModEq ∘ CartesianToKeplerian
    
    # Get inverse
    inv_composed = inv(composed)
    @test inv_composed isa ComposedTransformation
    
    # Apply forward and inverse
    modeq_coord = composed(cart_coord, μ)
    cart_back = inv_composed(modeq_coord, μ)
    
    @test params(cart_back) ≈ params(cart_coord) atol=1e-6 rtol=1e-8
end

@testset "compose() function" begin
    μ = 398600.4418
    cart_coord = Cartesian(7000.0, 1000.0, 2000.0, 1.0, 7.0, -0.5)
    
    # Test compose function (equivalent to ∘)
    composed = compose(KeplerianToModEq, CartesianToKeplerian)
    @test composed isa ComposedTransformation
    
    modeq_coord = composed(cart_coord, μ)
    @test modeq_coord isa ModEq
end

@testset "compose() with IdentityTransformation" begin
    μ = 398600.4418
    
    id_trans = IdentityTransformation()
    cart_to_kep = CartesianToKeplerian
    
    # compose(trans, IdentityTransformation) should return trans
    composed1 = compose(cart_to_kep, id_trans)
    @test composed1 === cart_to_kep
    
    # compose(IdentityTransformation, trans) should return trans
    composed2 = compose(id_trans, cart_to_kep)
    @test composed2 === cart_to_kep
    
    # compose(IdentityTransformation, IdentityTransformation) should return IdentityTransformation
    composed3 = compose(id_trans, id_trans)
    @test composed3 isa IdentityTransformation
end

@testset "Transformation inverse: Base.inv()" begin
    # Test that all transformation pairs have proper inverses
    
    # Cartesian ↔ Keplerian
    @test inv(CartesianToKeplerian) === KeplerianToCartesian
    @test inv(KeplerianToCartesian) === CartesianToKeplerian
    
    # Cartesian ↔ Spherical
    @test inv(CartesianToSpherical) === SphericalToCartesian
    @test inv(SphericalToCartesian) === CartesianToSpherical
    
    # Cartesian ↔ Cylindrical
    @test inv(CartesianToCylindrical) === CylindricalToCartesian
    @test inv(CylindricalToCartesian) === CartesianToCylindrical
    
    # Cartesian ↔ Delaunay
    @test inv(CartesianToDelaunay) === DelaunayToCartesian
    @test inv(DelaunayToCartesian) === CartesianToDelaunay
    
    # Cartesian ↔ Milankovich
    @test inv(CartesianToMilankovich) === MilankovichToCartesian
    @test inv(MilankovichToCartesian) === CartesianToMilankovich
    
    # Keplerian ↔ ModEq
    @test inv(KeplerianToModEq) === ModEqToKeplerian
    @test inv(ModEqToKeplerian) === KeplerianToModEq
    
    # Keplerian ↔ USM7
    @test inv(KeplerianToUSM7) === USM7ToKeplerian
    @test inv(USM7ToKeplerian) === KeplerianToUSM7
    
    # USM7 ↔ USM6
    @test inv(USM7ToUSM6) === USM6ToUSM7
    @test inv(USM6ToUSM7) === USM7ToUSM6
    
    # USM7 ↔ USMEM
    @test inv(USM7ToUSMEM) === USMEMToUSM7
    @test inv(USMEMToUSM7) === USM7ToUSMEM
end

@testset "Multiple composed transformations" begin
    μ = 398600.4418
    cart_coord = Cartesian(7000.0, 1000.0, 2000.0, 1.0, 7.0, -0.5)
    
    # Create a long chain: Cartesian → Keplerian → USM7 → USM6
    chain = USM7ToUSM6 ∘ KeplerianToUSM7 ∘ CartesianToKeplerian
    
    usm6_coord = chain(cart_coord, μ)
    @test usm6_coord isa USM6
    
    # Test inverse chain
    inv_chain = inv(chain)
    cart_back = inv_chain(usm6_coord, μ)
    @test params(cart_back) ≈ params(cart_coord) atol=1e-6 rtol=1e-8
end