using Test
using AstroCoords
using LinearAlgebra
using Random

# Load uncertainty propagation utilities (from examples)
include("../../examples/probabilistic_programming/uncertainty_propagation.jl")

@testset "Uncertainty Propagation Methods" begin
    μ = 398600.4418  # Earth gravitational parameter (km³/s²)

    @testset "Linear propagation" begin
        # Create Keplerian orbit with small uncertainty in semi-major axis
        kep = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.5)
        Σ_kep = Diagonal([100.0, 1e-6, 1e-8, 1e-8, 1e-8, 1e-8])  # 10 km std in a
        unc_kep = UncertainCoord(kep, CovarianceUncertainty(Σ_kep))

        # Transform to Cartesian
        transform = (k, μ) -> Cartesian(k, μ)
        unc_cart = propagate_linear(unc_kep, transform, μ)

        # Check that output is UncertainCoord with Cartesian
        @test unc_cart isa UncertainCoord
        @test unc_cart.coord isa Cartesian
        @test unc_cart.uncertainty isa CovarianceUncertainty

        # Check output covariance is positive semi-definite
        Σ_cart = unc_cart.uncertainty.cov
        @test issymmetric(Σ_cart)
        λ = eigvals(Σ_cart)
        @test all(λ .>= -1e-10)  # Eigenvalues should be non-negative (allowing small numerical error)

        # Check nominal transformation is correct
        cart_nominal = Cartesian(kep, μ)
        @test unc_cart.coord.x ≈ cart_nominal.x atol=1e-10
        @test unc_cart.coord.y ≈ cart_nominal.y atol=1e-10
        @test unc_cart.coord.z ≈ cart_nominal.z atol=1e-10
    end

    @testset "Linear propagation - round trip" begin
        # Keplerian -> Cartesian -> Keplerian should preserve covariance structure
        kep = Keplerian(7000.0, 0.1, 0.2, 0.3, 0.4, 0.5)
        Σ_kep = Diagonal([100.0, 0.001, 1e-4, 1e-4, 1e-4, 1e-4])
        unc_kep = UncertainCoord(kep, CovarianceUncertainty(Σ_kep))

        # Forward and back
        to_cart = (k, μ) -> Cartesian(k, μ)
        unc_cart = propagate_linear(unc_kep, to_cart, μ)

        to_kep = (c, μ) -> Keplerian(c, μ)
        unc_kep_back = propagate_linear(unc_cart, to_kep, μ)

        # Nominal coordinates should match (within tolerance for numerical round-trip)
        @test unc_kep_back.coord.a ≈ kep.a atol=1e-8 rtol=1e-10
        @test unc_kep_back.coord.e ≈ kep.e atol=1e-10 rtol=1e-10
        @test unc_kep_back.coord.i ≈ kep.i atol=1e-10 rtol=1e-10

        # Covariance should be similar (but not exact due to linearization errors)
        # Just check it's the right size and positive definite
        @test size(unc_kep_back.uncertainty.cov) == size(Σ_kep)
        @test issymmetric(unc_kep_back.uncertainty.cov)
        @test all(eigvals(unc_kep_back.uncertainty.cov) .>= -1e-10)
    end

    @testset "Unscented Transform - sigma point generation" begin
        mean = [1.0, 2.0, 3.0]
        cov = Diagonal([0.1, 0.2, 0.3])

        sigma_points, weights_mean, weights_cov = generate_sigma_points(mean, cov)

        # Should have 2n+1 points
        n = length(mean)
        @test size(sigma_points) == (n, 2n+1)

        # First point should be mean
        @test sigma_points[:, 1] ≈ mean atol=1e-12

        # Weights should sum to 1
        @test sum(weights_mean) ≈ 1.0 atol=1e-12
        @test sum(weights_cov) ≈ 1.0 atol=1e-12

        # Reconstruct mean and covariance from sigma points (should match input)
        mean_recon, cov_recon = reconstruct_statistics(
            sigma_points, weights_mean, weights_cov
        )
        @test mean_recon ≈ mean atol=1e-10
        @test cov_recon ≈ cov atol=1e-10
    end

    @testset "Unscented Transform propagation" begin
        kep = Keplerian(7000.0, 0.1, 0.2, 0.3, 0.4, 0.5)
        Σ_kep = Diagonal([100.0, 0.01, 1e-3, 1e-3, 1e-3, 1e-3])
        unc_kep = UncertainCoord(kep, CovarianceUncertainty(Σ_kep))

        transform = (k, μ) -> Cartesian(k, μ)
        unc_cart = propagate_unscented(unc_kep, transform, μ)

        # Check output structure
        @test unc_cart isa UncertainCoord
        @test unc_cart.coord isa Cartesian
        @test unc_cart.uncertainty isa CovarianceUncertainty

        # Check covariance properties
        Σ_cart = unc_cart.uncertainty.cov
        @test issymmetric(Σ_cart)
        @test all(eigvals(Σ_cart) .>= -1e-10)

        # Compare with linear propagation (should be similar for small uncertainties)
        unc_cart_linear = propagate_linear(unc_kep, transform, μ)

        # Means should be close
        @test unc_cart.coord.x ≈ unc_cart_linear.coord.x atol=1.0 rtol=0.01
        @test unc_cart.coord.y ≈ unc_cart_linear.coord.y atol=1.0 rtol=0.01
        @test unc_cart.coord.z ≈ unc_cart_linear.coord.z atol=1.0 rtol=0.01
    end

    @testset "Monte Carlo propagation" begin
        Random.seed!(12345)  # For reproducibility

        kep = Keplerian(7000.0, 0.05, 0.1, 0.2, 0.3, 0.5)
        Σ_kep = Diagonal([50.0, 0.005, 1e-3, 1e-3, 1e-3, 1e-3])
        unc_kep = UncertainCoord(kep, CovarianceUncertainty(Σ_kep))

        transform = (k, μ) -> Cartesian(k, μ)
        unc_cart, samples = propagate_monte_carlo(unc_kep, transform, μ; n_samples=10000)

        # Check output structure
        @test unc_cart isa UncertainCoord
        @test unc_cart.coord isa Cartesian
        @test size(samples) == (6, 10000)

        # Check covariance properties
        Σ_cart = unc_cart.uncertainty.cov
        @test issymmetric(Σ_cart)
        @test all(eigvals(Σ_cart) .>= -1e-10)

        # Compare with unscented transform (should be similar)
        unc_cart_ut = propagate_unscented(unc_kep, transform, μ)

        # Means should be close (within 3σ / √N ≈ 3% for N=10000)
        @test unc_cart.coord.x ≈ unc_cart_ut.coord.x atol=10.0 rtol=0.05
        @test unc_cart.coord.y ≈ unc_cart_ut.coord.y atol=10.0 rtol=0.05
        @test unc_cart.coord.z ≈ unc_cart_ut.coord.z atol=10.0 rtol=0.05

        # Diagonal covariance elements should be similar
        for i in 1:6
            @test Σ_cart[i, i] ≈ unc_cart_ut.uncertainty.cov[i, i] atol=Σ_cart[i, i]*0.1 rtol=0.1
        end
    end

    @testset "Edge case - zero uncertainty" begin
        kep = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.5)
        Σ_zero = zeros(6, 6)
        unc_kep = UncertainCoord(kep, CovarianceUncertainty(Σ_zero))

        transform = (k, μ) -> Cartesian(k, μ)

        # Linear propagation should work
        unc_cart = propagate_linear(unc_kep, transform, μ)
        @test all(abs.(unc_cart.uncertainty.cov) .< 1e-10)

        # Unscented should also work (degenerate case)
        # Note: This tests numerical stability when covariance is singular
        @test_throws Union{LinearAlgebra.PosDefException,ErrorException} propagate_unscented(
            unc_kep, transform, μ
        )
    end

    @testset "Edge case - large uncertainty" begin
        kep = Keplerian(7000.0, 0.5, 0.5, 0.5, 0.5, 0.5)
        # Very large uncertainty (unrealistic but tests numerical robustness)
        Σ_large = Diagonal([1e6, 0.1, 0.5, 0.5, 0.5, 0.5])
        unc_kep = UncertainCoord(kep, CovarianceUncertainty(Σ_large))

        transform = (k, μ) -> Cartesian(k, μ)

        # Should not error (even if results are not physically meaningful)
        unc_cart_linear = propagate_linear(unc_kep, transform, μ)
        @test unc_cart_linear isa UncertainCoord
        @test issymmetric(unc_cart_linear.uncertainty.cov)

        unc_cart_ut = propagate_unscented(unc_kep, transform, μ)
        @test unc_cart_ut isa UncertainCoord
        @test issymmetric(unc_cart_ut.uncertainty.cov)
    end
end
