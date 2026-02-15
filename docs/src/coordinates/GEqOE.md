# Generalized Equinoctial Orbital Elements

The Generalized Equinoctial Orbital Elements (GEqOE) extend the classical equinoctial elements to account for perturbations derived from a disturbing potential. Unlike the standard modified equinoctial elements, which are only valid under Keplerian dynamics, GEqOE absorbs the effect of conservative perturbations (e.g., gravitational harmonics) directly into the element definitions. This makes them particularly well-suited for high-fidelity orbit propagation and uncertainty quantification when the dominant perturbations are potential-based.

The formulation replaces the classical two-body integrals (energy, eccentricity vector, angular momentum) with their generalized counterparts that include the disturbing potential U, yielding elements that vary more slowly under perturbed motion.

## Components

The GEqOE state vector `[ν, p₁, p₂, L, q₁, q₂]` consists of six elements:

* **Orbital Shape and Size**
    * generalized mean motion (ν): Defined as ν = (1/μ)(-2E)^(3/2), where E is the total energy including the disturbing potential. This generalizes the Keplerian mean motion.
    * first generalized eccentricity component (p₁): p₁ = g sin(Ψ), where g is the magnitude of the generalized Laplace vector and Ψ is the generalized longitude of pericenter.
    * second generalized eccentricity component (p₂): p₂ = g cos(Ψ), analogous to the equinoctial eccentricity components but defined using the generalized eccentricity.

* **Position within the Orbit**
    * generalized mean longitude (L): Replaces the true longitude and encodes the satellite's position along its orbit in a generalized sense.

* **Orbital Orientation**
    * first inclination component (q₁): q₁ = tan(i/2) sin(Ω), identical to the classical equinoctial definition.
    * second inclination component (q₂): q₂ = tan(i/2) cos(Ω), identical to the classical equinoctial definition.

## Computed Properties

The following derived quantities are available on a `GEqOE` instance:

* `g`: Generalized eccentricity magnitude, √(p₁² + p₂²)
* `Ψ`: Generalized longitude of pericenter, atan(p₁, p₂)

Additional helper functions:

* `L₀(geq, t)`: Generalized mean longitude at epoch, L₀ = L - νt
* `computed_a(geq, μ)`: Generalized semi-major axis, a = (μ/ν²)^(1/3)
* `computed_ρ(geq, μ)`: Generalized semi-latus rectum, ρ = a(1 - g²)

## Extension Requirement

Unlike most coordinate sets in AstroCoords.jl, GEqOE transformations require a dynamics model to evaluate the disturbing potential at the current state. This means the `AstroForceModels.jl` extension must be loaded to perform coordinate transformations:

```julia
using AstroCoords
using AstroForceModels
using ComponentArrays

μ = 3.986004415e5
cart_state = Cartesian(-1076.225, -6765.896, -332.309, 9.357, -3.312, -1.188)

# Keplerian dynamics (U = 0, equivalent to modified equinoctial)
dynamics = CentralBodyDynamicsModel(KeplerianGravityAstroModel(μ=μ))
p = ComponentVector(JD=2.460310e6)
t = 0.0

geqoe_state = GEqOE(cart_state, μ, dynamics, p, t)

# Round-trip back to Cartesian
cart_roundtrip = Cartesian(geqoe_state, μ, dynamics, p, t)
```

For non-Keplerian dynamics with gravitational harmonics:

```julia
using SatelliteToolboxGravityModels
using SatelliteToolboxTransformations

eop_data = fetch_iers_eop()
grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

dynamics = CentralBodyDynamicsModel(GravityHarmonicsAstroModel(;
    gravity_model=grav_coeffs, eop_data=eop_data, order=4, degree=4
))

geqoe_state = GEqOE(cart_state, μ, dynamics, p, t)
```

## Singularities

* Non-singular for circular orbits (g = 0)
* Non-singular for equatorial orbits (i = 0)
* Singular for retrograde equatorial orbits (i = π) and rectilinear motion
* Requires negative total energy (E < 0), i.e., bound orbits only

## References
[1]: Baù, G., Hernando-Ayuso, J., & Bombardelli, C. (2021). "A generalization of the equinoctial orbital elements." Celestial Mechanics and Dynamical Astronomy, 133(9), 1-32.
