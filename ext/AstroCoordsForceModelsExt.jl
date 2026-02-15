module AstroCoordsForceModelsExt

using AstroCoords
using AstroForceModels: AstroForceModels
using AstroForceModels: AbstractDynamicsModel
using ComponentArrays
using StaticArrays
using LinearAlgebra

# Keplerian reference model used to isolate the disturbing potential
const _KeplerianGravityModel = AstroForceModels.KeplerianGravityAstroModel

"""
    _disturbing_potential(u, μ, model, p, t)

Compute the disturbing potential U (paper Eq. 3) by subtracting the Keplerian
central-body potential from the total gravitational potential returned by the
force model:  U = V_total - V_keplerian.
"""
@inline function _disturbing_potential(u, μ, model, p, t)
    V_total = AstroForceModels.potential(u, p, t, model.gravity_model)
    V_kepler = AstroForceModels.potential(u, p, t, _KeplerianGravityModel(μ=μ))
    return V_total - V_kepler
end

# ~~~~~~~~~~~~~~~ Generalized Kepler Solver ~~~~~~~~~~~~~~~ #

"""
    _solve_kepler_generalized(L, p₁, p₂; tol=1e-14, maxiter=50)

Solve the generalized Kepler equation (Eq. 25):
    L = K + p₁ cos(K) - p₂ sin(K)

Uses Newton-Raphson iteration.
"""
function _solve_kepler_generalized(
    L::Number,
    p₁::Number,
    p₂::Number;
    tol::Float64=1e-14,
    maxiter::Int=200,
    verbosity::Val{verbose}=Val(false),
) where {verbose}
    K = L
    for _ in 1:maxiter
        sK, cK = sincos(K)
        f = L - K - p₁ * cK + p₂ * sK
        fp = -1 + p₁ * sK + p₂ * cK
        δK = -f / fp
        K += δK
        abs(δK) < tol && return K
    end
    verbosity &&
        @warn "Generalized Kepler equation did not converge after $maxiter iterations"
    return K
end

# ~~~~~~~~~~~~~~~ cart2geqoe ~~~~~~~~~~~~~~~ #

"""
    cart2geqoe(u, μ, model, p, t)

Transform Cartesian coordinates to Generalized Equinoctial Orbital Elements (GEqOE).

The disturbing potential U is extracted from the central body gravity model
(paper Eq. 3: F = P - ∇U).

# Arguments
- `u::AbstractVector{<:Number}`: Cartesian state vector [x; y; z; ẋ; ẏ; ż]
- `μ::Number`: Gravitational parameter
- `model::AbstractDynamicsModel`: Dynamics model (e.g., `CentralBodyDynamicsModel`)
- `p::ComponentVector`: Simulation parameters
- `t::Number`: Current simulation time [s]

# Returns
- `SVector{6}`: GEqOE state vector [ν; p₁; p₂; L; q₁; q₂]
"""
function AstroCoords.cart2geqoe(
    u::AbstractVector{T}, μ::V, model::AbstractDynamicsModel, p::ComponentVector, t::TT
) where {T<:Number,V<:Number,TT<:Number}
    # Disturbing potential U (paper Eq. 3): total gravity minus Keplerian central body
    U = _disturbing_potential(u, μ, model, p, t)
    W = typeof(U)
    RT = promote_type(T, V, W, TT)

    # Extract position and velocity
    r = SVector{3,RT}(u[1], u[2], u[3])
    ṙ = SVector{3,RT}(u[4], u[5], u[6])

    r_mag = norm(r)
    ṙ_dot = dot(r, ṙ) / r_mag  # Radial velocity

    # Angular momentum vector and magnitude
    h_vec = cross(r, ṙ)
    h = norm(h_vec)

    # Effective potential energy: Ueff = h²/(2r²) + U (Section 2.1)
    U_eff = h^2 / (2 * r_mag^2) + U

    # Generalized angular momentum: c = √(2r²Ueff) (Eq. 4)
    c = sqrt(2 * r_mag^2 * U_eff)

    # Total energy: E = ½|v̇|² - μ/r + U (Section 2.1)
    E = RT(0.5) * dot(ṙ, ṙ) - μ / r_mag + U

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
    geqoe2cart(u, μ, model, p, t)

Transform GEqOE to Cartesian coordinates using an AstroForceModels dynamics model.

The angular momentum h is recovered from the generalized angular momentum c via
h = √(c² - 2r²U), where U is the central body gravity potential evaluated at the
reconstructed position (paper Section 4, Eq. 41).

# Arguments
- `u::AbstractVector{<:Number}`: GEqOE state vector [ν; p₁; p₂; L; q₁; q₂]
- `μ::Number`: Gravitational parameter
- `model::AbstractDynamicsModel`: Dynamics model (e.g., `CentralBodyDynamicsModel`)
- `p::ComponentVector`: Simulation parameters
- `t::Number`: Current simulation time [s]

# Returns
- `SVector{6}`: Cartesian state vector [x; y; z; ẋ; ẏ; ż]
"""
function AstroCoords.geqoe2cart(
    u::AbstractVector{T}, μ::V, model::AbstractDynamicsModel, p::ComponentVector, t::TT
) where {T<:Number,V<:Number,TT<:Number}
    RT = promote_type(T, V, TT)

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

    # Evaluate disturbing potential U at reconstructed position (Section 4)
    # h = √(c² - 2r²U), where U is the disturbing (non-Keplerian) potential
    u_full = SVector{6}(r[1], r[2], r[3], zero(RT), zero(RT), zero(RT))
    U = _disturbing_potential(u_full, μ, model, p, t)
    h = sqrt(c^2 - 2 * r_mag^2 * U)

    # Velocity vector: v = ṙ_dot·eᵣ + (h/r)·eₐ (Eq. 41)
    v = ṙ_dot * eᵣ + (h / r_mag) * eₐ

    return SVector{6,RT}(r[1], r[2], r[3], v[1], v[2], v[3])
end

# ~~~~~~~~~~~~~~~ Transformation Callable Methods ~~~~~~~~~~~~~~~ #
# These provide the actual implementations for the structs defined in transformations.jl

function (::AstroCoords.CartesianToGEqOETransform)(
    x::AstroCoords.Cartesian,
    μ::Number,
    model::AbstractDynamicsModel,
    p::ComponentVector,
    t::Number,
)
    geqoe_vec = AstroCoords.cart2geqoe(AstroCoords.params(x), μ, model, p, t)
    return AstroCoords.GEqOE(geqoe_vec...)
end

function (::AstroCoords.GEqOEToCartesianTransform)(
    x::AstroCoords.GEqOE,
    μ::Number,
    model::AbstractDynamicsModel,
    p::ComponentVector,
    t::Number,
)
    cart_vec = AstroCoords.geqoe2cart(AstroCoords.params(x), μ, model, p, t)
    return AstroCoords.Cartesian(cart_vec...)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
# Direct: Cartesian <-> GEqOE (more specific than auto-generated args... versions)

function AstroCoords.GEqOE(
    x::AstroCoords.Cartesian{T},
    μ::Number,
    model::AbstractDynamicsModel,
    p::ComponentVector,
    t::Number,
) where {T<:Number}
    geqoe_vec = AstroCoords.cart2geqoe(AstroCoords.params(x), μ, model, p, t)
    return AstroCoords.GEqOE(geqoe_vec...)
end

function AstroCoords.Cartesian(
    x::AstroCoords.GEqOE{T},
    μ::Number,
    model::AbstractDynamicsModel,
    p::ComponentVector,
    t::Number,
) where {T<:Number}
    cart_vec = AstroCoords.geqoe2cart(AstroCoords.params(x), μ, model, p, t)
    return AstroCoords.Cartesian(cart_vec...)
end

end # module
