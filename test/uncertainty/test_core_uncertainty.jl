using Test
using AstroCoords
using LinearAlgebra
using StaticArrays

@testset "Core Uncertainty Types" begin
    @testset "CovarianceUncertainty construction" begin
        # Valid covariance matrix
        Σ = Diagonal([1.0, 2.0, 3.0])
        unc = CovarianceUncertainty(Σ)
        @test unc isa CovarianceUncertainty{Float64}
        @test unc.cov ≈ Σ atol=1e-15

        # Asymmetric matrix should error
        Σ_asym = [1.0 2.0; 3.0 4.0]
        @test_throws ArgumentError CovarianceUncertainty(Σ_asym)

        # Symmetric matrix (non-diagonal)
        Σ_sym = [1.0 0.5; 0.5 2.0]
        unc_sym = CovarianceUncertainty(Σ_sym)
        @test unc_sym isa CovarianceUncertainty{Float64}
        @test issymmetric(unc_sym.cov)
    end

    @testset "UncertainCoord construction" begin
        # Keplerian with covariance uncertainty
        kep = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.5)
        Σ = Diagonal([100.0, 0.001, 1e-6, 1e-6, 1e-6, 1e-6])
        unc_kep = UncertainCoord(kep, CovarianceUncertainty(Σ))

        @test unc_kep isa UncertainCoord
        @test unc_kep.coord === kep
        @test unc_kep.uncertainty isa CovarianceUncertainty
        @test unc_kep.uncertainty.cov ≈ Σ atol=1e-15

        # Cartesian with uncertainty
        cart = Cartesian(6524.834, 6862.875, 6448.296, 4.901327, 5.533756, -1.976341)
        Σ_cart = Diagonal([1.0, 1.0, 1.0, 0.01, 0.01, 0.01])  # km, km/s
        unc_cart = UncertainCoord(cart, CovarianceUncertainty(Σ_cart))

        @test unc_cart.coord === cart
        @test unc_cart.uncertainty.cov ≈ Σ_cart atol=1e-15
    end

    @testset "UncertainCoord property access" begin
        kep = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.5)
        Σ = Diagonal([100.0, 0.001, 1e-6, 1e-6, 1e-6, 1e-6])
        unc_kep = UncertainCoord(kep, CovarianceUncertainty(Σ))

        # Should access coordinate properties
        @test unc_kep.a ≈ 7000.0 atol=1e-12
        @test unc_kep.e ≈ 0.01 atol=1e-12
        @test unc_kep.i ≈ 0.1 atol=1e-12

        # Should access uncertainty
        @test unc_kep.uncertainty isa CovarianceUncertainty
    end

    @testset "Type stability" begin
        kep_f64 = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.5)
        Σ_f64 = Diagonal([100.0, 0.001, 1e-6, 1e-6, 1e-6, 1e-6])
        unc_f64 = UncertainCoord(kep_f64, CovarianceUncertainty(Σ_f64))

        @test @inferred(UncertainCoord(kep_f64, CovarianceUncertainty(Σ_f64))) isa
            UncertainCoord

        kep_f32 = Keplerian(7000.0f0, 0.01f0, 0.1f0, 0.2f0, 0.3f0, 0.5f0)
        Σ_f32 = Diagonal(Float32[100.0, 0.001, 1e-6, 1e-6, 1e-6, 1e-6])
        unc_f32 = UncertainCoord(kep_f32, CovarianceUncertainty(Σ_f32))

        @test eltype(unc_f32.coord) === Float32
        @test eltype(unc_f32.uncertainty.cov) === Float32
    end

    @testset "Edge cases" begin
        # Zero covariance (perfect knowledge)
        kep = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.5)
        Σ_zero = zeros(6, 6)
        unc_zero = UncertainCoord(kep, CovarianceUncertainty(Σ_zero))
        @test unc_zero.uncertainty.cov ≈ Σ_zero atol=1e-15

        # Large uncertainty
        Σ_large = 1e6 * Diagonal(ones(6))
        unc_large = UncertainCoord(kep, CovarianceUncertainty(Σ_large))
        @test unc_large.uncertainty.cov ≈ Σ_large atol=1e-9

        # Non-diagonal covariance (correlated uncertainties)
        Σ_corr = [
            100.0 10.0 0.0 0.0 0.0 0.0;
            10.0 0.001 0.0 0.0 0.0 0.0;
            0.0 0.0 1e-6 0.0 0.0 0.0;
            0.0 0.0 0.0 1e-6 0.0 0.0;
            0.0 0.0 0.0 0.0 1e-6 0.0;
            0.0 0.0 0.0 0.0 0.0 1e-6
        ]
        unc_corr = UncertainCoord(kep, CovarianceUncertainty(Σ_corr))
        @test issymmetric(unc_corr.uncertainty.cov)
        @test unc_corr.uncertainty.cov[1, 2] ≈ 10.0 atol=1e-15
    end
end
