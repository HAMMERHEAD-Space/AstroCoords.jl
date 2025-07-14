export cart2StiefelScheifele, StiefelScheifele2cart, set_stiefelscheifele_configurations

"""
    StiefelScheifele2cart(ss_vec, μ; ϕ₀, t₀, DU, TU, W, flag_time)

Converts a Stiefel-Scheifele state vector to a Cartesian state vector.

# Arguments
- `ss_vec::AbstractVector`: The Stiefel-Scheifele state vector.
- `μ::Number`: The gravitational parameter of the central body.
- `ϕ₀::Number`: The initial fictitious time.
- `t₀::Number`: The initial physical time.
- `DU::Number`: The distance unit.
- `TU::Number`: The time unit.
- `W::Number`: The perturbing potential energy.
- `flag_time::AbstractTimeType`: The time element representation (defaults to `PhysicalTime()`).

# Returns
- `Cartesian`: The Cartesian state vector.
"""
function StiefelScheifele2cart(
    ss_vec::AbstractVector,
    μ::Number;
    ϕ₀::Number,
    t₀::Number,
    DU::Number,
    TU::Number,
    W::Number,
    flag_time::AbstractTimeType=PhysicalTime(),
)
    α = SVector{4}(ss_vec[1], ss_vec[2], ss_vec[3], ss_vec[4])
    β = SVector{4}(ss_vec[5], ss_vec[6], ss_vec[7], ss_vec[8])
    ω = ss_vec[9]

    ##################################################
    #* 1. Auxiliary Quantities
    ##################################################
    sph2, cph2 = sincos(0.5 * ϕ₀)

    ##################################################
    #* 2. Position in the Inertial Frame
    ##################################################
    #SS State Vector and its derivative
    u = α * cph2 + β * sph2
    du = 0.5 * (-α * sph2 + β * cph2)

    # Position in inertial frame
    r = SVector{3}(
        u[1]^2 - u[2]^2 - u[3]^2 + u[4]^2,
        2 * (u[1] * u[2] - u[3] * u[4]),
        2 * (u[1] * u[3] + u[2] * u[4]),
    )
    rmag = u[1]^2 + u[2]^2 + u[3]^2 + u[4]^2

    ##################################################
    #* 3. Velocity in the Inertial Frame
    ##################################################
    v = SVector{3}(
        (u[1] * du[1] - u[2] * du[2] - u[3] * du[3] + u[4] * du[4]) * (4 * ω / rmag),
        (u[2] * du[1] + u[1] * du[2] - u[4] * du[3] - u[3] * du[4]) * (4 * ω / rmag),
        (u[3] * du[1] + u[4] * du[2] + u[1] * du[3] + u[2] * du[4]) * (4 * ω / rmag),
    )

    return SVector{6}(
        r[1] * DU, r[2] * DU, r[3] * DU, v[1] * DU / TU, v[2] * DU / TU, v[3] * DU / TU
    )
end

"""
    cart2StiefelScheifele(cart_vec, μ; DU, TU, W, ϕ, t₀, flag_time)

Converts Cartesian coordinates to Stiefel-Scheifele coordinates.

# Arguments
- `cart_vec::AbstractVector`: The Cartesian state vector.
- `μ::Number`: The gravitational parameter of the central body.
- `DU::Number`: The distance unit.
- `TU::Number`: The time unit.
- `W::Number`: The perturbing potential energy.
- `ϕ₀::Number`: The initial fictitious time.
- `t₀::Number`: The initial physical time.
- `flag_time::AbstractTimeType`: The time element representation (defaults to `PhysicalTime()`).

# Returns
- `StiefelScheifele`: The Stiefel-Scheifele state vector.
"""
function cart2StiefelScheifele(
    cart_vec::AbstractVector,
    μ::Number;
    DU::Number,
    TU::Number,
    W::Number,
    ϕ₀::Number,
    t₀::Number,
    flag_time::AbstractTimeType=PhysicalTime(),
)
    ##################################################
    #* 1. Non-Dimensionalize
    ##################################################
    x = SVector{3}(cart_vec[1] / DU, cart_vec[2] / DU, cart_vec[3] / DU)
    xdot = SVector{3}(
        cart_vec[4] / (DU / TU), cart_vec[5] / (DU / TU), cart_vec[6] / (DU / TU)
    )
    pot0 = W / (DU^2 / TU^2)
    Ksq = μ / (DU^3 / TU^2)

    ##################################################
    #* 2. Auxiliary Quantities
    ##################################################
    Rmag = norm(x)
    sph2, cph2 = sincos(ϕ₀ / 2)

    ##################################################
    #* 2. Compute Components
    ##################################################
    totEn0 = 0.5 * dot(xdot, xdot) - Ksq / Rmag + pot0
    ω = √(-totEn0 / 2)

    if x[1] >= 0
        u₁ = 0.0
        u₄ = √(0.5 * (Rmag + x[1]))
        u₂ = (x[2]*u₁ + x[3]*u₄) / (Rmag + x[1])
        u₃ = (x[3]*u₁ - x[2]*u₄) / (Rmag + x[1])
    else
        u₂ = 0.0
        u₃ = √(0.5 * (Rmag - x[1]))
        u₁ = (x[2]*u₂ + x[3]*u₃) / (Rmag - x[1])
        u₄ = (x[3]*u₂ - x[2]*u₃) / (Rmag - x[1])
    end

    # Derivatives of the K-S parameters wrt the independent variable

    du₁ = (u₁ * xdot[1] + u₂ * xdot[2] + u₃ * xdot[3]) / (4 * ω)
    du₂ = (-u₂ * xdot[1] + u₁ * xdot[2] + u₄ * xdot[3]) / (4 * ω)
    du₃ = (-u₃ * xdot[1] - u₄ * xdot[2] + u₁ * xdot[3]) / (4 * ω)
    du₄ = (u₄ * xdot[1] - u₃ * xdot[2] + u₂ * xdot[3]) / (4 * ω)

    α = cph2 * SVector{4}(u₁, u₂, u₃, u₄) - 2 * sph2 * SVector{4}(du₁, du₂, du₃, du₄)
    β = sph2 * SVector{4}(u₁, u₂, u₃, u₄) + 2 * cph2 * SVector{4}(du₁, du₂, du₃, du₄)

    # Time / Time Element
    if flag_time isa PhysicalTime
        t_val = t₀ * TU
    elseif flag_time isa LinearTime
        t_val =
            t₀ * TU + dot(SVector{4}(u₁, u₂, u₃, u₄), SVector{4}(du₁, du₂, du₃, du₄)) / ω
    end

    return SVector{10}(α[1], α[2], α[3], α[4], β[1], β[2], β[3], β[4], ω, t_val)
end

"""
    set_stiefelscheifele_configurations(state::Cartesian, μ; ϕ=0.0, W=0.0, t₀=0.0, flag_time=0)

Computes the necessary parameters for Stiefel-Scheifele transformations based on an
initial Cartesian state.

# Arguments
- `state::Cartesian`: The initial Cartesian state vector (position and velocity).
- `μ`: The gravitational parameter of the central body.
- `ϕ`: (Optional) The initial fictitious time (defaults to 0.0).
- `W`: (Optional) The perturbing potential energy (defaults to 0.0).
- `t₀`: (Optional) The initial physical time (defaults to 0.0).
- `flag_time`: (Optional) Flag to determine time element representation (0 for physical time, 1 for linear time element; defaults to 0).

# Returns
- A `NamedTuple` containing `DU`, `TU`, `W`, `ϕ`, `t₀`, and `flag_time`.
"""
function set_stiefelscheifele_configurations(
    state::AbstractVector,
    μ::Number;
    W::Number=0.0,
    t₀::Number=0.0,
    flag_time::AbstractTimeType=PhysicalTime(),
)
    DU = norm(SVector{3}(state[1], state[2], state[3]))
    TU = sqrt(DU^3 / μ)
    ϕ₀ = computeϕ₀(state, μ, DU, TU, W)
    return (; DU, TU, W, ϕ₀, t₀, flag_time)
end

"""
    get_stiefelscheifele_time(u::AbstractVector{T}, flag_time::AbstractTimeType) where {T<:Number}

Computes the time element for Stiefel-Scheifele transformations based on an
initial Cartesian state.

# Arguments
- `u::AbstractVector{T}`: The Stiefel-Scheifele state vector.

# Keyword Arguments
- `t₀::Number`: The initial physical time.
- `TU::Number`: The time unit.
- `flag_time::AbstractTimeType`: The time element representation (defaults to `PhysicalTime()`).

# Returns
- The time element.
"""
function get_stiefelscheifele_time(
    u::AbstractVector; ϕ::Number, flag_time::AbstractTimeType=PhysicalTime()
)
    if flag_time isa PhysicalTime
        return u[10]
    elseif flag_time isa LinearTime
        sph, cph = sincos(ϕ)
        αsq = u[1]^2 + u[2]^2 + u[3]^2 + u[4]^2
        βsq = u[5]^2 + u[6]^2 + u[7]^2 + u[8]^2
        αβ = u[1] * u[5] + u[2] * u[6] + u[3] * u[7] + u[4] * u[8]
        return u[10] + 0.5 * ((αsq - βsq) / 2.0 * sph - αβ * cph) / u[9]
    end
end
