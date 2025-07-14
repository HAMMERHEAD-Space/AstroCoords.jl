# Stiefel-Scheifele

The Stiefel-Scheifele (SS) transformation is a method used to regularize the equations of motion in the two-body problem, similar to the Kustaanheimo-Stiefel (KS) transformation. It is particularly effective for orbits with high eccentricity, where the standard Newtonian formulation can suffer from numerical instability near the pericenter. By transforming the problem into a higher-dimensional space with a fictitious time, the SS formulation yields a set of linear, constant-coefficient differential equations, which are more stable to integrate numerically.

## Components

The Stiefel-Scheifele state vector is composed of ten elements that describe the regularized orbit:

*   **α Vector (α₁ - α₄)**: A four-component vector representing the transformed position. It is derived from a linear combination of the KS position vector and its derivative, effectively encoding the orbital geometry in a regularized form.

*   **β Vector (β₁ - β₄)**: A four-component vector that is also a linear combination of the KS position and velocity vectors. Together with the α vector, it defines the state of the system in the transformed phase space.

*   **Total Energy Parameter (ω)**: A parameter directly related to the total energy of the orbit, where `ω = sqrt(-E/2)`. This formulation ensures that the parameter remains real and positive for bound orbits.

*   **Time Element (t)**: Represents time within the formulation, which can be either physical time or a linear time element, depending on the specific implementation and the needs of the integration scheme.

## References
[1]: Stiefel, E. L. and Scheifele, G. "Linear and Regular Celestial Mechanics", Springer-Verlag, 1971. 
[2]: Amato, Davide. "THALASSA: Orbit propagator for near-Earth and cislunar space." Astrophysics Source Code Library (2019): ascl-1905. 