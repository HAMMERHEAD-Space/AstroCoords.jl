# Poincaré

Poincaré canonical orbital elements are a set of canonical variables widely used in Hamiltonian celestial mechanics and perturbation theory. They are derived from Delaunay action-angle variables via a canonical transformation that replaces the eccentricity and inclination action-angle pairs with Cartesian-style coordinate-momentum pairs. This transformation eliminates the singularities that arise in classical orbital elements for circular orbits (e → 0) and equatorial orbits (i → 0).

## Components

The Poincaré elements consist of three canonical conjugate pairs:

* **Semi-Major Axis Pair**
    * Canonical Action (Λ): Related to the semi-major axis by Λ = √(μ|a|), where μ is the gravitational parameter.
    * Mean Longitude (λ): The sum of the mean anomaly, argument of periapsis, and RAAN: λ = M + ω + Ω.

* **Eccentricity Pair (Cartesian-style)**
    * ξ = √(2P) cos(ω̃): Eccentricity cosine component, where ω̃ = ω + Ω is the longitude of periapsis.
    * η = -√(2P) sin(ω̃): Eccentricity sine component.
    * The eccentricity action P depends on orbit type:
        * Elliptic (a > 0): P = Λ - G = Λ(1 - √(1 - e²)) ∈ [0, Λ]
        * Hyperbolic (a < 0): P = Λ + G = Λ(1 + √(e² - 1)) ∈ (Λ, ∞)

* **Inclination Pair (Cartesian-style)**
    * p = √(2Q) cos(Ω): Inclination cosine component, where Q = G(1 - cos(i)) is the inclination action and G is the Delaunay angular momentum.
    * q = -√(2Q) sin(Ω): Inclination sine component.

## Relation to Delaunay Variables

The Poincaré elements are obtained from the Delaunay variables (L, G, H, l, g, h) through the canonical transformation:

| Poincaré | Delaunay |
|----------|----------|
| Λ = L | L = √(μ\|a\|) |
| λ = l + g + h | Mean longitude |
| P = L - G (elliptic) or L + G (hyperbolic) | Eccentricity action |
| ω̃ = -(g + h) | Negative longitude of periapsis |
| Q = G - H | Inclination action |
| Ω = -h | Negative RAAN |

The Cartesian-style pairs (ξ, η) and (p, q) are then defined as:
- (ξ, η) = √(2P) × (cos(-ω̃), sin(-ω̃))
- (p, q) = √(2Q) × (cos(-Ω), sin(-Ω))

The orbit type is recovered unambiguously during inversion by comparing P against Λ:
- P ≤ Λ → elliptic (G = Λ - P, a = Λ²/μ)
- P > Λ → hyperbolic (G = P - Λ, a = -Λ²/μ)

## Properties

* **Non-singular** for circular orbits (e → 0): As e → 0, both ξ and η smoothly approach zero regardless of the undefined longitude of periapsis.
* **Non-singular** for equatorial orbits (i → 0): As i → 0, both p and q smoothly approach zero regardless of the undefined RAAN.
* **Canonical**: The transformation preserves the symplectic structure, making these elements ideal for Hamiltonian perturbation theory.

## References
[1]: Murray, C.D. and Dermott, S.F. "Solar System Dynamics." Cambridge University Press (1999).
[2]: Laskar, J. and Robutel, P. "Stability of the Planetary Three-Body Problem." Celestial Mechanics and Dynamical Astronomy 62 (1995): 193-217.
[3]: Morbidelli, A. "Modern Celestial Mechanics: Aspects of Solar System Dynamics." Taylor & Francis (2002).
