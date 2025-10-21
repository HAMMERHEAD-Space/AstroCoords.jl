# EDromo

The EDromo state vector is a non-singular set of 8 orbital elements designed for accurate trajectory propagation, especially in the presence of perturbations. It combines aspects of equinoctial and quaternion elements to avoid singularities associated with near-circular or near-equatorial orbits. The formulation is valid for all orbit types (elliptic, parabolic, and hyperbolic) and uses regularization techniques that make propagation exceptionally numerically stable.

## Components

The EDromo state vector `[ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈]` consists of eight elements:

*   **In-Plane Elements**
    *   `ζ₁`: In-plane element related to eccentricity and argument of periapsis.
    *   `ζ₂`: In-plane element related to eccentricity and argument of periapsis.
    *   `ζ₃`: Related to the orbital energy. Its sign defines the type of orbit:
        *   `ζ₃ > 0` for elliptic orbits (negative energy).
        *   `ζ₃ = 0` for parabolic orbits (zero energy).
        *   `ζ₃ < 0` for hyperbolic orbits (positive energy).

*   **Quaternion Elements (Orientation)**
    *   `ζ₄`: Quaternion element 1.
    *   `ζ₅`: Quaternion element 2.
    *   `ζ₆`: Quaternion element 3.
    *   `ζ₇`: Quaternion element 4.

*   **Time Element**
    *   `ζ₈`: Time element, which can represent physical, constant, or linear time.

## Special Configurations

Unlike other coordinate systems, the transformation to and from `EDromo` coordinates requires a set of special configuration parameters. These are necessary for non-dimensionalization and for defining the fictitious time frame.

A helper function is provided to generate these parameters from a standard Cartesian state.

## References
[1]: Baù, G., Bombardelli, C., Peláez, J., and Lorenzini, E., "Nonsingular orbital elements for special perturbations in the two-body problem". MNRAS 454(3), pp. 2890-2908. 2015. 
[2]: Amato, Davide. "THALASSA: Orbit propagator for near-Earth and cislunar space." Astrophysics Source Code Library (2019): ascl-1905.