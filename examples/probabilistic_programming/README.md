# Probabilistic Programming Examples for AstroCoords.jl

This directory contains examples demonstrating AstroCoords.jl's compatibility with probabilistic programming frameworks for uncertainty propagation through coordinate transformations.

## Setup

These examples require additional dependencies not included in the main AstroCoords.jl package. To use these examples:

```julia
using Pkg
Pkg.activate("examples/probabilistic_programming")
Pkg.instantiate()
```

## Examples

### 1. Linear Uncertainty Propagation (`01_linear_propagation.jl`)

Demonstrates first-order covariance propagation using Jacobian matrices. This is the standard method used in spacecraft navigation.

**Reference**: Tapley, B. D., Schutz, B. E., & Born, G. H. (2004). *Statistical Orbit Determination*. Academic Press. Chapter 4.

**Method**: Σ_out = J * Σ_in * J^T, where J is the Jacobian of the transformation.

### 2. Unscented Transform (`02_unscented_transform.jl`)

Implements the Unscented Transform for nonlinear uncertainty propagation using sigma points.

**Reference**: Julier, S. J., & Uhlmann, J. K. (1997). "New extension of the Kalman filter to nonlinear systems." *Signal Processing, Sensor Fusion, and Target Recognition VI*, 3068, 182-193.

**Method**: Propagate 2n+1 sigma points through nonlinear transformation, then reconstruct statistics.

### 3. Monte Carlo with Gaussian Mixture Models (`03_monte_carlo_gmm.jl`)

Shows how to handle non-Gaussian uncertainties using Gaussian Mixture Models and Monte Carlo sampling.

**Method**: Represent uncertainty as weighted sum of Gaussians, sample and propagate through transformations.

### 4. Turing.jl Integration (`04_turing_example.jl`)

Demonstrates using Turing.jl for Bayesian inference with orbital mechanics. Shows how to:
- Define priors over orbital elements
- Condition on observations
- Sample posterior distributions

### 5. RxInfer.jl Integration (`05_rxinfer_example.jl`)

Shows belief propagation on factor graphs for coordinate transformations. Useful for:
- Real-time state estimation
- Sensor fusion
- Online learning

### 6. Stan Integration (`06_stan_example.jl`)

Demonstrates calling Julia functions from Stan models for hybrid inference.

## Core Library Addition

A minimal `UncertainCoord` type has been added to the core library in `src/uncertainty.jl`:

```julia
using AstroCoords
using LinearAlgebra

# Keplerian orbit with uncertainty
kep = Keplerian(7000.0, 0.01, 0.1, 0.2, 0.3, 0.5)
Σ = Diagonal([100.0, 0.001, 0.0, 0.0, 0.0, 0.0])  # uncertainty in a and e
unc_kep = UncertainCoord(kep, CovarianceUncertainty(Σ))
```

For advanced uncertainty propagation methods, use the examples in this directory.

## References

1. Tapley, B. D., Schutz, B. E., & Born, G. H. (2004). *Statistical Orbit Determination*. Academic Press.
2. Julier, S. J., & Uhlmann, J. K. (1997). "New extension of the Kalman filter to nonlinear systems." *Signal Processing, Sensor Fusion, and Target Recognition VI*, 3068, 182-193.
3. Carpenter, J., Clifford, P., & Fearnhead, P. (1999). "Improved particle filter for nonlinear problems." *IEE Proceedings-Radar, Sonar and Navigation*, 146(1), 2-7.