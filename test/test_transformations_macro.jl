using Test
using AstroCoords

@testset "Transformation Macro and Generated Functions" begin
    @testset "Macro Expansion - Transformation Pairs" begin
        # Test that macro-generated types exist
        # The macro should generate GCRFToITRF and ITRFToGCRF types
        @test isdefined(AstroCoords, :GCRFToITRF)
        @test isdefined(AstroCoords, :ITRFToGCRF)
    end

    @testset "Forward Transformation Constructor" begin
        # Test GCRFToITRF constructor
        t = 0.0  # Some time value
        
        trans = AstroCoords.GCRFToITRF(t)
        
        @test trans isa AstroCoords.GCRFToITRF
        @test typeof(trans.time) == Float64
    end

    @testset "Reverse Transformation Constructor" begin
        # Test ITRFToGCRF constructor
        t = 0.0
        
        trans = AstroCoords.ITRFToGCRF(t)
        
        @test trans isa AstroCoords.ITRFToGCRF
        @test typeof(trans.time) == Float64
    end

    @testset "Transformation with Timestamp" begin
        # Test that transformations work with time-dependent paths
        using Dates
        
        # Create transformation at specific time
        t = DateTime(2024, 1, 1, 12, 0, 0)
        trans_gcrf_to_itrf = AstroCoords.GCRFToITRF(t)
        
        @test trans_gcrf_to_itrf.time == t
    end

    @testset "Round-Trip Transformation" begin
        # Test that forward then reverse gets back to original
        # (identity up to numerical precision)
        t = 0.0
        
        # Create test state in GCRF
        gcrf_state = AstroCoords.Cartesian(
            7000.0, 0.0, 0.0,  # position (km)
            0.0, 7.5, 0.0      # velocity (km/s)
        )
        
        # Forward: GCRF → ITRF
        trans_forward = AstroCoords.GCRFToITRF(t)
        itrf_state = trans_forward(gcrf_state)
        
        # Reverse: ITRF → GCRF
        trans_reverse = AstroCoords.ITRFToGCRF(t)
        gcrf_back = trans_reverse(itrf_state)
        
        # Check round-trip preservation
        @test gcrf_back.x ≈ gcrf_state.x atol=1e-10
        @test gcrf_back.y ≈ gcrf_state.y atol=1e-10
        @test gcrf_back.z ≈ gcrf_state.z atol=1e-10
        @test gcrf_back.ẋ ≈ gcrf_state.ẋ atol=1e-10
        @test gcrf_back.ẏ ≈ gcrf_state.ẏ atol=1e-10
        @test gcrf_back.ż ≈ gcrf_state.ż atol=1e-10
    end

    @testset "Transformation Type Stability" begin
        # Test that transformations are type-stable
        t = 0.0
        trans = AstroCoords.GCRFToITRF(t)
        state = AstroCoords.Cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0)
        
        # Type inference should work
        @inferred trans(state)
    end

    @testset "Transformation Different Time Types" begin
        # Test with different time representations
        using Dates
        
        # Float time
        trans1 = AstroCoords.GCRFToITRF(0.0)
        @test trans1 isa AstroCoords.GCRFToITRF
        
        # DateTime
        trans2 = AstroCoords.GCRFToITRF(DateTime(2024, 1, 1))
        @test trans2 isa AstroCoords.GCRFToITRF
        
        # Integer time
        trans3 = AstroCoords.GCRFToITRF(0)
        @test trans3 isa AstroCoords.GCRFToITRF
    end

    @testset "Multiple Transformation Pairs" begin
        # If other frame pairs are defined, test them
        # This tests the macro's ability to generate multiple pairs
        t = 0.0
        
        # Test GCRF-ITRF pair exists
        @test isdefined(AstroCoords, :GCRFToITRF)
        @test isdefined(AstroCoords, :ITRFToGCRF)
        
        # Can create instances
        @test AstroCoords.GCRFToITRF(t) isa AstroCoords.GCRFToITRF
        @test AstroCoords.ITRFToGCRF(t) isa AstroCoords.ITRFToGCRF
    end
end
