using Test
using AstroCoords

@testset "Module Loading and Initialization" begin
    @testset "Module Loads Successfully" begin
        # Test that the module loads without errors
        @test isdefined(AstroCoords, :Cartesian)
        @test isdefined(AstroCoords, :Keplerian)
    end
    
    @testset "All Coordinate Types Exported" begin
        # Verify main coordinate types are accessible
        @test isdefined(AstroCoords, :Cartesian)
        @test isdefined(AstroCoords, :Keplerian)
        @test isdefined(AstroCoords, :Cylindrical)
        @test isdefined(AstroCoords, :Spherical)
        @test isdefined(AstroCoords, :Delaunay)
        @test isdefined(AstroCoords, :Milankovich)
        @test isdefined(AstroCoords, :ModEq)
        @test isdefined(AstroCoords, :J2EqOE)
        @test isdefined(AstroCoords, :USM6)
        @test isdefined(AstroCoords, :USM7)
        @test isdefined(AstroCoords, :USMEM)
        @test isdefined(AstroCoords, :EDromo)
        @test isdefined(AstroCoords, :KustaanheimoStiefel)
        @test isdefined(AstroCoords, :StiefelScheifele)
    end
    
    @testset "Transformation Functions Exported" begin
        # Verify transformation functions are accessible
        @test isdefined(AstroCoords, :CartesianToKeplerian)
        @test isdefined(AstroCoords, :KeplerianToCartesian)
        @test isdefined(AstroCoords, :CartesianToModEq)
        @test isdefined(AstroCoords, :ModEqToCartesian)
    end
    
    @testset "Utility Functions Exported" begin
        # Verify utility functions are accessible
        @test isdefined(AstroCoords, :angle_between_vectors)
        @test isdefined(AstroCoords, :meanMotion)
        @test isdefined(AstroCoords, :orbitalPeriod)
        @test isdefined(AstroCoords, :KeplerSolver)
    end
    
    @testset "COORDINATE_SET_ALIASES Dict Populated" begin
        # Test that the aliases dictionary exists and is populated
        @test isdefined(AstroCoords, :COORDINATE_SET_ALIASES)
        aliases = AstroCoords.COORDINATE_SET_ALIASES
        
        @test haskey(aliases, "Cartesian")
        @test haskey(aliases, "Keplerian")
        @test haskey(aliases, "Cylindrical")
        @test haskey(aliases, "Spherical")
        @test haskey(aliases, "Delaunay")
        @test haskey(aliases, "EDromo")
        @test haskey(aliases, "J2EqOE")
        @test haskey(aliases, "KustaanheimoStiefel")
        @test haskey(aliases, "StiefelScheifele")
        @test haskey(aliases, "Milankovich")
        @test haskey(aliases, "ModifiedEquinoctial")
        @test haskey(aliases, "USM6")
        @test haskey(aliases, "USM7")
        @test haskey(aliases, "USMEM")
        
        # Verify aliases map to correct types
        @test aliases["Cartesian"] === Cartesian
        @test aliases["Keplerian"] === Keplerian
        @test aliases["ModifiedEquinoctial"] === ModEq
    end
    
    @testset "No Namespace Collisions" begin
        # Test that exported names don't collide with Base/Core
        # (This is a safety check to ensure module doesn't break existing code)
        @test AstroCoords.params !== Base.params  # Should be module-specific
        @test AstroCoords.meanMotion !== Base.mean  # Different function
    end
    
    @testset "Core Types Defined" begin
        # Verify abstract types are defined
        @test isdefined(AstroCoords, :AstroCoord)
        @test isdefined(AstroCoords, :Transformation)
        @test isdefined(AstroCoords, :AstrodynamicsTransformation)
    end
    
    @testset "Anomaly Conversions Exported" begin
        # Verify anomaly conversion functions
        @test isdefined(AstroCoords, :trueAnomaly2MeanAnomaly)
        @test isdefined(AstroCoords, :meanAnomaly2TrueAnomaly)
        @test isdefined(AstroCoords, :trueAnomaly2EccentricAnomaly)
        @test isdefined(AstroCoords, :eccentricAnomaly2TrueAnomaly)
    end
    
    @testset "Regularized Coordinate Support" begin
        # Verify regularized coordinate configuration is available
        @test isdefined(AstroCoords, :RegularizedCoordinateConfig)
        @test isdefined(AstroCoords, :PhysicalTime)
        @test isdefined(AstroCoords, :ConstantTime)
        @test isdefined(AstroCoords, :LinearTime)
    end
end
