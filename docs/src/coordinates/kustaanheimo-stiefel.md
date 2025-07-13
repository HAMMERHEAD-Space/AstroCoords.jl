# Kustaanheimo-Stiefel

The Kustaanheimo-Stiefel (KS) transformation regularizes the two-body problem by transforming the 3D Cartesian coordinates into a 4D space, which eliminates the singularity at `r=0`. This is particularly useful for modeling highly elliptical or near-rectilinear orbits where the orbiting body passes very close to the primary. The KS state vector also includes the negative of the total energy and a time element, which can be physical or fictitious time.

## Components

The Kustaanheimo-Stiefel state vector consists of ten elements:

*   **KS Position Vector (u₁ - u₄)**: Four components that represent the transformed position in a 4D space. This transformation is what regularizes the equations of motion.

*   **KS Velocity Vector (u₁' - u₄')**: The four corresponding velocity components in the transformed 4D space. These are the derivatives of the KS position components with respect to a fictitious time.

*   **Negative Total Energy (h)**: The negative of the total orbital energy (`-E`). This value is conserved in the two-body problem and is used in the KS equations.

*   **Time Element (τ)**: A time-like variable that can represent either physical time (`t`) or a linear time element, depending on the chosen formulation. This provides flexibility in the integration scheme.

## References
[1]: Stiefel, E. L. and Scheifele, G. "Linear and Regular Celestial Mechanics", Springer-Verlag, 1971. 
[2]: Amato, Davide. "THALASSA: Orbit propagator for near-Earth and cislunar space." Astrophysics Source Code Library (2019): ascl-1905.