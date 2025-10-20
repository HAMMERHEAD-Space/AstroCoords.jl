# Delaunay

The Delaunay Orbital Elements are an alternative formulation of classical orbital elements expressed in terms of action-angle variables. These elements are widely used in perturbation theory and celestial mechanics because they simplify the equations of motion by transforming them into canonical Hamiltonian form. This makes them particularly useful for studying long-term orbital behavior and complex gravitational interactions.

![image](../assets/delaunay.jpg)
*Image of Delaunay Orbital Elements [1]*

## Components

The Delaunay variables are called action-angle coordinates as the Hamiltonian only depends on the action variables.

* **Action Coordinates**
    * Canonical Keplerian Energy (L): Instaneous energy of the Keplerian orbit.
    * Canonical Total Angular Momentum (G): Instaneous total angular momentum of the Keplerian orbit.
    * Canonical Normal Angular Momentum (Relative to Equator) (H): Instaneous normal angular momentum of the Keplerian orbit.
* **Angle Coordinates**
    * Mean Anomaly (M): The angle between the satellite's current position and the periapsis, measured from the central body at a specific moment in time treated as thought the orbit was circular.
    * Argument of periapsis (ω): The angle between the ascending node and the point of closest approach to the central body (periapsis). It defines the orientation of the orbit within the plane.
    * Longitude of the ascending node (Ω): The angle from a fixed reference direction (typically the vernal equinox) to the point where the satellite crosses the equatorial plane from south to north (ascending node).

## Exotic Properties

In addition to the six standard Delaunay elements stored as fields, the `Delaunay` coordinate type provides computed properties for additional anomaly types:

* **Eccentric Anomaly (E)**: Computed property representing an auxiliary angle used in orbital mechanics calculations. Automatically calculated from the mean anomaly and eccentricity (derived from the L and G action variables).
* **True Anomaly (f)**: Computed property that represents the actual angular position of the satellite relative to periapsis. Also calculated from the mean anomaly and derived eccentricity.

These exotic properties are accessible as regular properties (e.g., `del.E`, `del.f`) but are computed on-demand rather than stored. The eccentricity required for these calculations is automatically extracted from the relationship e = √(1 - (G/L)²), ensuring consistency with the canonical action variables while providing convenient access to commonly needed derived quantities.

## References
[1]: Laskar, Jacques. "Andoyer construction for Hill and Delaunay variables." Celestial Mechanics and Dynamical Astronomy 128.4 (2017): 475-482.
[2]: https://www.researchgate.net/publication/2143998_The_Averaged_Dynamics_of_the_Hydrogen_Atom_in_Crossed_Electric_and_Magnetic_Fields_as_a_Perturbed_Kepler_Problem
