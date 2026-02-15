# ~~~~~~~~~~~~~~~ Generalized Kepler Solver ~~~~~~~~~~~~~~~ #

"""
    _solve_kepler_generalized(L, p₁, p₂; tol=1e-14, maxiter=200)

Solve the generalized Kepler equation (Eq. 25):
    L = K + p₁ cos(K) - p₂ sin(K)

Uses Newton-Raphson iteration.
"""
function _solve_kepler_generalized(
    L::Number, p₁::Number, p₂::Number; tol::Float64=1e-14, maxiter::Int=200
)
    K = L
    for _ in 1:maxiter
        sK, cK = sincos(K)
        f = L - K - p₁ * cK + p₂ * sK
        fp = -1 + p₁ * sK + p₂ * cK
        δK = -f / fp
        K += δK
        abs(δK) < tol && return K
    end
    return K
end

# ~~~~~~~~~~~~~~~ cart2geqoe ~~~~~~~~~~~~~~~ #

"""
    cart2geqoe(u, μ, config::RegularizedCoordinateConfig)

Transform Cartesian coordinates to Generalized Equinoctial Orbital Elements (GEqOE).

The perturbing potential W is taken from the `config` struct. For a pure Keplerian
orbit, set `W = 0`. For perturbed orbits, precompute the disturbing potential
(V_total - V_keplerian) externally and pass it via the config.

Reference: Baù, G., Hernando-Ayuso, J., & Bombardelli, C. (2021).
"A generalization of the equinoctial orbital elements."
Celestial Mechanics and Dynamical Astronomy, 133(9), 1-32.

# Arguments
- `u::AbstractVector{<:Number}`: Cartesian state vector [x; y; z; ẋ; ẏ; ż]
- `μ::Number`: Gravitational parameter
- `config::RegularizedCoordinateConfig`: Configuration containing the perturbing potential `W`

# Returns
- `SVector{6}`: GEqOE state vector [ν; p₁; p₂; L; q₁; q₂]
"""
function cart2geqoe(
    u::AbstractVector{T}, μ::V, config::RegularizedCoordinateConfig
) where {T<:Number,V<:Number}
    W = config.W
    RT = promote_type(T, V, typeof(W))

    # Extract position and velocity
    r = SVector{3,RT}(u[1], u[2], u[3])
    ṙ = SVector{3,RT}(u[4], u[5], u[6])

    r_mag = norm(r)
    ṙ_dot = dot(r, ṙ) / r_mag  # Radial velocity

    # Angular momentum vector and magnitude
    h_vec = cross(r, ṙ)
    h = norm(h_vec)

    # Effective potential energy: Ueff = h²/(2r²) + W (Section 2.1)
    U_eff = h^2 / (2 * r_mag^2) + W

    # Generalized angular momentum: c = √(2r²Ueff) (Eq. 4)
    c = sqrt(2 * r_mag^2 * U_eff)

    # Total energy: E = ½|v̇|² - μ/r + W (Section 2.1)
    E = RT(0.5) * dot(ṙ, ṙ) - μ / r_mag + W

    # Element 1: Generalized mean motion ν = (1/μ)(-2E)^(3/2) (Eq. 12)
    if E >= 0
        error("Total energy E must be negative for GEqOE (E = $E)")
    end
    ν = (1 / μ) * (-2 * E)^(RT(3) / 2)

    # Generalized semi-major axis: a = -μ/(2E) (Eq. 6)
    a = -μ / (2 * E)

    # Generalized semi-latus rectum: ρ = c²/μ (from Eqs. 15-16)
    ρ = c^2 / μ

    # Elements 5 & 6: q₁, q₂ (Eqs. 30-31)
    q₁ = h_vec[1] / (h + h_vec[3])
    q₂ = -h_vec[2] / (h + h_vec[3])

    # Equinoctial frame unit vectors (Eq. 33)
    γ = 1 + q₁^2 + q₂^2
    eₓ = SVector{3,RT}((1 - q₁^2 + q₂^2) / γ, (2 * q₁ * q₂) / γ, (-2 * q₁) / γ)
    eᵧ = SVector{3,RT}((2 * q₁ * q₂) / γ, (1 + q₁^2 - q₂^2) / γ, (2 * q₂) / γ)

    # True longitude
    eᵣ = r / r_mag
    cosL = dot(eᵣ, eₓ)
    sinL = dot(eᵣ, eᵧ)

    # Elements 2 & 3: p₁, p₂ (Eqs. 34-35)
    p₁ = (ρ / r_mag - 1) * sinL - (c * ṙ_dot / μ) * cosL
    p₂ = (ρ / r_mag - 1) * cosL + (c * ṙ_dot / μ) * sinL

    # Generalized eccentric longitude K (Eq. 36, Appendix A)
    w = sqrt(μ / a)
    S = (μ + c * w - r_mag * ṙ_dot^2) * sinL - ṙ_dot * (c + w * r_mag) * cosL
    C = (μ + c * w - r_mag * ṙ_dot^2) * cosL + ṙ_dot * (c + w * r_mag) * sinL
    K = atan(S, C)

    # Element 4: Generalized mean longitude (Eq. 25, compact form)
    L = K + (C * p₁ - S * p₂) / (μ + c * w)

    return SVector{6,RT}(ν, p₁, p₂, L, q₁, q₂)
end

# ~~~~~~~~~~~~~~~ geqoe2cart ~~~~~~~~~~~~~~~ #

"""
    geqoe2cart(u, μ, config::RegularizedCoordinateConfig)

Transform GEqOE to Cartesian coordinates.

The perturbing potential W is taken from the `config` struct. It must be the same
value used in the forward transformation.

# Arguments
- `u::AbstractVector{<:Number}`: GEqOE state vector [ν; p₁; p₂; L; q₁; q₂]
- `μ::Number`: Gravitational parameter
- `config::RegularizedCoordinateConfig`: Configuration containing the perturbing potential `W`

# Returns
- `SVector{6}`: Cartesian state vector [x; y; z; ẋ; ẏ; ż]
"""
function geqoe2cart(
    u::AbstractVector{T}, μ::V, config::RegularizedCoordinateConfig
) where {T<:Number,V<:Number}
    W = config.W
    RT = promote_type(T, V, typeof(W))

    # Extract GEqOE elements
    ν, p₁, p₂, L, q₁, q₂ = u

    # Derived quantities (Eqs. 15-17)
    a = (μ / ν^2)^(RT(1) / 3)
    g² = p₁^2 + p₂^2
    ρ = a * (1 - g²)
    c = sqrt(μ * ρ)

    # Solve generalized Kepler equation for K (Eq. 25)
    K = _solve_kepler_generalized(L, p₁, p₂)
    sinK, cosK = sincos(K)

    # Orbital distance and radial velocity (Eqs. 26-27)
    r_mag = a * (1 - p₁ * sinK - p₂ * cosK)
    ṙ_dot = sqrt(μ * a) / r_mag * (p₂ * sinK - p₁ * cosK)

    # True longitude from eccentric longitude (Eq. 38)
    α = 1 / (1 + sqrt(1 - g²))
    sinL = (a / r_mag) * (α * p₁ * p₂ * cosK + (1 - α * p₂^2) * sinK - p₁)
    cosL = (a / r_mag) * (α * p₁ * p₂ * sinK + (1 - α * p₁^2) * cosK - p₂)

    # Equinoctial frame unit vectors (Eq. 33)
    γ = 1 + q₁^2 + q₂^2
    eₓ = SVector{3,RT}((1 - q₁^2 + q₂^2) / γ, (2 * q₁ * q₂) / γ, (-2 * q₁) / γ)
    eᵧ = SVector{3,RT}((2 * q₁ * q₂) / γ, (1 + q₁^2 - q₂^2) / γ, (2 * q₂) / γ)

    # Orbital frame unit vectors (Eq. 40)
    eᵣ = eₓ * cosL + eᵧ * sinL
    eₐ = eᵧ * cosL - eₓ * sinL

    # Position vector: r = r·eᵣ (Eq. 41)
    r = r_mag * eᵣ

    # Angular momentum from generalized angular momentum (Section 4)
    # h = √(c² - 2r²W)
    h = sqrt(c^2 - 2 * r_mag^2 * W)

    # Velocity vector: v = ṙ_dot·eᵣ + (h/r)·eₐ (Eq. 41)
    v = ṙ_dot * eᵣ + (h / r_mag) * eₐ

    return SVector{6,RT}(r[1], r[2], r[3], v[1], v[2], v[3])
end
