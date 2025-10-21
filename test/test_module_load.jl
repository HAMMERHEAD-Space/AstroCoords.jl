using Test
using AstroCoords

@testset "Module Load and Initialization" begin
    @testset "Module loads successfully" begin
        # Test that the module itself is defined
        @test isdefined(Main, :AstroCoords)
        @test AstroCoords isa Module
    end

    @testset "COORDINATE_SET_ALIASES Dictionary" begin
        # Test that COORDINATE_SET_ALIASES is defined and populated
        @test isdefined(AstroCoords, :COORDINATE_SET_ALIASES)
        aliases = AstroCoords.COORDINATE_SET_ALIASES
        
        # Test that it's a dictionary
        @test aliases isa Dict
        
        # Test that it's not empty
        @test !isempty(aliases)
        
        # Test specific expected entries
        @test haskey(aliases, "Cartesian")
        @test haskey(aliases, "Keplerian")
        @test haskey(aliases, "ModifiedEquinoctial")
        @test haskey(aliases, "EDromo")
        @test haskey(aliases, "KustaanheimoStiefel")
        @test haskey(aliases, "StiefelScheifele")
        @test haskey(aliases, "Cylindrical")
        @test haskey(aliases, "Spherical")
        @test haskey(aliases, "Delaunay")
        @test haskey(aliases, "Milankovich")
        @test haskey(aliases, "J2EqOE")
        @test haskey(aliases, "USM6")
        @test haskey(aliases, "USM7")
        @test haskey(aliases, "USMEM")
        
        # Test that aliases map to correct types
        @test aliases["Cartesian"] === Cartesian
        @test aliases["Keplerian"] === Keplerian
        @test aliases["ModifiedEquinoctial"] === ModEq
        @test aliases["EDromo"] === EDromo
        @test aliases["KustaanheimoStiefel"] === KustaanheimoStiefel
        @test aliases["StiefelScheifele"] === StiefelScheifele
        @test aliases["Cylindrical"] === Cylindrical
        @test aliases["Spherical"] === Spherical
        @test aliases["Delaunay"] === Delaunay
        @test aliases["Milankovich"] === Milankovich
        @test aliases["J2EqOE"] === J2EqOE
        @test aliases["USM6"] === USM6
        @test aliases["USM7"] === USM7
        @test aliases["USMEM"] === USMEM
        
        # Test that all mapped types are actually coordinate types
        for (name, coord_type) in aliases
            @test coord_type <: AstroCoord
        end
    end

    @testset "Core Exports Available" begin
        # Test that key types and functions are exported and available
        @test isdefined(AstroCoords, :Cartesian)
        @test isdefined(AstroCoords, :Keplerian)
        @test isdefined(AstroCoords, :ModEq)
        @test isdefined(AstroCoords, :EDromo)
        @test isdefined(AstroCoords, :KustaanheimoStiefel)
        @test isdefined(AstroCoords, :StiefelScheifele)
        
        # Test transformation types
        @test isdefined(AstroCoords, :CartesianToKeplerian)
        @test isdefined(AstroCoords, :KeplerianToCartesian)
        
        # Test utility functions
        @test isdefined(AstroCoords, :params)
        @test isdefined(AstroCoords, :meanMotion)
        @test isdefined(AstroCoords, :orbitalPeriod)
        
        # Test config types
        @test isdefined(AstroCoords, :RegularizedCoordinateConfig)
        @test isdefined(AstroCoords, :PhysicalTime)
        @test isdefined(AstroCoords, :ConstantTime)
        @test isdefined(AstroCoords, :LinearTime)
    end

    @testset "Dependencies Loaded" begin
        # Test that required dependencies are available
        @test isdefined(AstroCoords, :LinearAlgebra)
        @test isdefined(AstroCoords, :StaticArrays)
        
        # Test that we can use types from dependencies
        using StaticArrays
        vec = SVector{3}(1.0, 2.0, 3.0)
        @test vec isa SVector
    end
end
