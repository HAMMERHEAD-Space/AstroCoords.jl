using LinearAlgebra
using StaticArrays

export meanMotion,
    orbitalPeriod, orbitalNRG, angularMomentumVector, angularMomentumQuantity

"""
    meanMotion(state::AstroCoord, μ)

Compute the mean motion of the orbit.

# Arguments
- `state::AstroCoord`: The orbital state.
- `μ`: Gravitational parameter.

# Returns
- Mean motion in rad/TU.
"""
function meanMotion(state::AstroCoord, μ)
    cart_state = Cartesian(state, μ)
    return meanMotion(cart_state, μ)
end

function meanMotion(state::Cartesian, μ)
    r = position(state)
    v = velocity(state)
    r_mag = norm(r)
    v_mag = norm(v)
    E = v_mag^2 / 2 - μ / r_mag
    a = -μ / (2 * E)
    return √(μ / a^3)
end

function meanMotion(state::Keplerian, μ)
    return √(μ / state.a^3)
end

# Specialized methods for regularized coordinates
function meanMotion(
    state::EDromo, μ, ϕ::Number, config::RegularizedCoordinateConfig
)
    # Convert to Cartesian using provided ϕ and config
    cart_state = Cartesian(state, μ, ϕ, config)
    return meanMotion(cart_state, μ)
end

function meanMotion(
    state::StiefelScheifele, μ, ϕ::Number, config::RegularizedCoordinateConfig
)
    # Convert to Cartesian using provided ϕ and config
    cart_state = Cartesian(state, μ, ϕ, config)
    return meanMotion(cart_state, μ)
end

"""
    orbitalPeriod(state::AstroCoord, μ)

Compute the orbital period.

# Arguments
- `state::AstroCoord`: The orbital state.
- `μ`: Gravitational parameter.

# Returns
- Orbital period in TU.
"""
function orbitalPeriod(state::AstroCoord, μ)
    return 2 * π / meanMotion(state, μ)
end

function orbitalPeriod(state::Keplerian, μ)
    return 2 * π * √(state.a^3 / μ)
end

# Specialized methods for regularized coordinates
function orbitalPeriod(
    state::EDromo, μ, ϕ::Number, config::RegularizedCoordinateConfig
)
    return 2 * π / meanMotion(state, μ, ϕ, config)
end

function orbitalPeriod(
    state::StiefelScheifele, μ, ϕ::Number, config::RegularizedCoordinateConfig
)
    return 2 * π / meanMotion(state, μ, ϕ, config)
end

"""
    orbitalNRG(state::AstroCoord, μ)

Compute the specific orbital energy.

# Arguments
- `state::AstroCoord`: The orbital state.
- `μ`: Gravitational parameter.

# Returns
- Specific orbital energy in (L²/T²).
"""
function orbitalNRG(state::AstroCoord, μ)
    cart_state = Cartesian(state, μ)
    return orbitalNRG(cart_state, μ)
end

function orbitalNRG(state::Cartesian, μ)
    r = position(state)
    v = velocity(state)
    r_mag = norm(r)
    v_mag = norm(v)
    return v_mag^2 / 2 - μ / r_mag
end

function orbitalNRG(state::Keplerian, μ)
    return -μ / (2 * state.a)
end

# Specialized methods for regularized coordinates
function orbitalNRG(
    state::EDromo, μ, ϕ::Number, config::RegularizedCoordinateConfig
)
    cart_state = Cartesian(state, μ, ϕ, config)
    return orbitalNRG(cart_state, μ)
end

function orbitalNRG(
    state::StiefelScheifele, μ, ϕ::Number, config::RegularizedCoordinateConfig
)
    cart_state = Cartesian(state, μ, ϕ, config)
    return orbitalNRG(cart_state, μ)
end

function orbitalNRG(
    state::KustaanheimoStiefel, μ, config::RegularizedCoordinateConfig
)
    cart_state = Cartesian(state, μ, config)
    return orbitalNRG(cart_state, μ)
end

"""
    angularMomentumVector(state::AstroCoord, μ)

Compute the angular momentum vector.

# Arguments
- `state::AstroCoord`: The orbital state.
- `μ`: Gravitational parameter.

# Returns
- Angular momentum vector as `SVector{3}`.
"""
function angularMomentumVector(state::AstroCoord, μ)
    cart_state = Cartesian(state, μ)
    return angularMomentumVector(cart_state, μ)
end

function angularMomentumVector(state::Cartesian, μ)
    r = position(state)
    v = velocity(state)
    return cross(r, v)
end

# Specialized methods for regularized coordinates
function angularMomentumVector(
    state::EDromo, μ, ϕ::Number, config::RegularizedCoordinateConfig
)
    cart_state = Cartesian(state, μ, ϕ, config)
    return angularMomentumVector(cart_state, μ)
end

function angularMomentumVector(
    state::StiefelScheifele, μ, ϕ::Number, config::RegularizedCoordinateConfig
)
    cart_state = Cartesian(state, μ, ϕ, config)
    return angularMomentumVector(cart_state, μ)
end

function angularMomentumVector(
    state::KustaanheimoStiefel, μ, config::RegularizedCoordinateConfig
)
    cart_state = Cartesian(state, μ, config)
    return angularMomentumVector(cart_state, μ)
end

"""
    angularMomentumQuantity(state::AstroCoord, μ)

Compute the magnitude of the angular momentum.

# Arguments
- `state::AstroCoord`: The orbital state.
- `μ`: Gravitational parameter.

# Returns
- Magnitude of angular momentum.
"""
function angularMomentumQuantity(state::AstroCoord, μ)
    return norm(angularMomentumVector(state, μ))
end

# Specialized methods for regularized coordinates
function angularMomentumQuantity(
    state::EDromo, μ, ϕ::Number, config::RegularizedCoordinateConfig
)
    return norm(angularMomentumVector(state, μ, ϕ, config))
end

function angularMomentumQuantity(
    state::StiefelScheifele, μ, ϕ::Number, config::RegularizedCoordinateConfig
)
    return norm(angularMomentumVector(state, μ, ϕ, config))
end

function angularMomentumQuantity(
    state::KustaanheimoStiefel, μ, config::RegularizedCoordinateConfig
)
    return norm(angularMomentumVector(state, μ, config))
end