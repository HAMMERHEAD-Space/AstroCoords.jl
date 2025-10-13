using Test
using AstroCoords
using StaticArrays

@testset "Cartesian Constructors" begin
    @testset "Inner Constructor" begin
        # Test basic inner constructor
        x, y, z = 7000.0, 0.0, 0.0
        ẋ, ẏ, ż = 0.0, 7.5, 0.0
        
        cart = AstroCoords.Cartesian{Float64}(x, y, z, ẋ, ẏ, ż)
        
        @test cart.x == x
        @test cart.y == y
        @test cart.z == z
        @test cart.ẋ == ẋ
        @test cart.ẏ == ẏ
        @test cart.ż == ż
    end

    @testset "Outer Constructor - Individual Values" begin
        # Test outer constructor with individual arguments
        x, y, z = 7000.0, 1000.0, 500.0
        ẋ, ẏ, ż = 1.0, 7.5, 0.5
        
        cart = AstroCoords.Cartesian(x, y, z, ẋ, ẏ, ż)
        
        @test cart.x == x
        @test cart.y == y
        @test cart.z == z
        @test cart.ẋ == ẋ
        @test cart.ẏ == ẏ
        @test cart.ż == ż
        @test typeof(cart) == AstroCoords.Cartesian{Float64}
    end

    @testset "Vector Constructor" begin
        # Test constructor from vector
        X = [7000.0, 1000.0, 500.0, 1.0, 7.5, 0.5]
        
        cart = AstroCoords.Cartesian(X)
        
        @test cart.x == X[1]
        @test cart.y == X[2]
        @test cart.z == X[3]
        @test cart.ẋ == X[4]
        @test cart.ẏ == X[5]
        @test cart.ż == X[6]
    end

    @testset "StaticVector Constructor - SVector" begin
        # Test with SVector
        X = SVector{6}(7000.0, 1000.0, 500.0, 1.0, 7.5, 0.5)
        
        cart = AstroCoords.Cartesian(X)
        
        @test cart.x == X[1]
        @test cart.y == X[2]
        @test cart.z == X[3]
        @test cart.ẋ == X[4]
        @test cart.ẏ == X[5]
        @test cart.ż == X[6]
    end

    @testset "StaticVector Constructor - MVector" begin
        # Test with MVector (mutable static array)
        X = MVector{6}(7000.0, 1000.0, 500.0, 1.0, 7.5, 0.5)
        
        cart = AstroCoords.Cartesian(X)
        
        @test cart.x == X[1]
        @test cart.y == X[2]
        @test cart.z == X[3]
        @test cart.ẋ == X[4]
        @test cart.ẏ == X[5]
        @test cart.ż == X[6]
    end

    @testset "Type Promotion - Float32 + Float64" begin
        # Test type promotion with mixed types
        x = Float32(7000.0)
        y = Float64(1000.0)
        z = Float32(500.0)
        ẋ = Float64(1.0)
        ẏ = Float32(7.5)
        ż = Float64(0.5)
        
        cart = AstroCoords.Cartesian(x, y, z, ẋ, ẏ, ż)
        
        # Should promote to Float64
        @test typeof(cart) == AstroCoords.Cartesian{Float64}
        @test cart.x ≈ 7000.0
        @test cart.y ≈ 1000.0
    end

    @testset "Type Promotion - Int + Float" begin
        # Test type promotion with integers and floats
        x = 7000
        y = 1000.0
        z = 500
        ẋ = 1
        ẏ = 7.5
        ż = 0
        
        cart = AstroCoords.Cartesian(x, y, z, ẋ, ẏ, ż)
        
        # Should promote to Float64
        @test typeof(cart) == AstroCoords.Cartesian{Float64}
        @test cart.x == 7000.0
        @test cart.ẏ == 7.5
    end

    @testset "Copy Constructor" begin
        # Test copy constructor
        cart1 = AstroCoords.Cartesian{Float64}(7000.0, 1000.0, 500.0, 1.0, 7.5, 0.5)
        cart2 = AstroCoords.Cartesian{Float64}(cart1)
        
        @test cart2.x == cart1.x
        @test cart2.y == cart1.y
        @test cart2.z == cart1.z
        @test cart2.ẋ == cart1.ẋ
        @test cart2.ẏ == cart1.ẏ
        @test cart2.ż == cart1.ż
    end

    @testset "Accessors" begin
        # Test field access
        cart = AstroCoords.Cartesian(7000.0, 1000.0, 500.0, 1.0, 7.5, 0.5)
        
        @test cart.x == 7000.0
        @test cart.y == 1000.0
        @test cart.z == 500.0
        @test cart.ẋ == 1.0
        @test cart.ẏ == 7.5
        @test cart.ż == 0.5
    end

    @testset "Type Stability" begin
        # Test type stability with @inferred
        x, y, z = 7000.0, 1000.0, 500.0
        ẋ, ẏ, ż = 1.0, 7.5, 0.5
        
        @inferred AstroCoords.Cartesian(x, y, z, ẋ, ẏ, ż)
        
        X = [x, y, z, ẋ, ẏ, ż]
        @inferred AstroCoords.Cartesian(X)
    end

    @testset "Equality" begin
        # Test that equivalent constructors produce equal objects
        x, y, z = 7000.0, 1000.0, 500.0
        ẋ, ẏ, ż = 1.0, 7.5, 0.5
        
        cart1 = AstroCoords.Cartesian(x, y, z, ẋ, ẏ, ż)
        cart2 = AstroCoords.Cartesian([x, y, z, ẋ, ẏ, ż])
        cart3 = AstroCoords.Cartesian(SVector{6}(x, y, z, ẋ, ẏ, ż))
        
        @test cart1.x == cart2.x == cart3.x
        @test cart1.y == cart2.y == cart3.y
        @test cart1.z == cart2.z == cart3.z
        @test cart1.ẋ == cart2.ẋ == cart3.ẋ
        @test cart1.ẏ == cart2.ẏ == cart3.ẏ
        @test cart1.ż == cart2.ż == cart3.ż
    end

    @testset "Zero State" begin
        # Test constructing zero state
        cart = AstroCoords.Cartesian(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        
        @test cart.x == 0.0
        @test cart.y == 0.0
        @test cart.z == 0.0
        @test cart.ẋ == 0.0
        @test cart.ẏ == 0.0
        @test cart.ż == 0.0
    end
end
