export get_stiefelscheifele_time

"""
    StiefelScheifele2cart(ss_vec, μ, ϕ::Number, config::RegularizedCoordinateConfig)

Converts a Stiefel-Scheifele state vector to a Cartesian state vector.

# Arguments
- `ss_vec::AbstractVector`: The Stiefel-Scheifele state vector.
- `μ::Number`: The gravitational parameter of the central body.
- `ϕ::Number`: Fictitious time parameter.
- `config::RegularizedCoordinateConfig`: Configuration parameters.

# Returns
- `Cartesian`: The Cartesian state vector.
"""
function StiefelScheifele2cart(
    ss_vec::AbstractVector, μ::Number, ϕ::Number, config::RegularizedCoordinateConfig
)
    DU, TU = config.DU, config.TU

    α = SVector{4}(ss_vec[1], ss_vec[2], ss_vec[3], ss_vec[4])
    β = SVector{4}(ss_vec[5], ss_vec[6], ss_vec[7], ss_vec[8])
    ω = ss_vec[9]

    ##################################################
    #* 1. Auxiliary Quantities
    ##################################################
    sph2, cph2 = sincos(0.5 * ϕ)

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
    cart2StiefelScheifele(cart_vec, μ, ϕ::Number, config::RegularizedCoordinateConfig)

Converts Cartesian coordinates to Stiefel-Scheifele coordinates.

# Arguments
- `cart_vec::AbstractVector`: The Cartesian state vector.
- `μ::Number`: The gravitational parameter of the central body.
- `ϕ::Number`: Fictitious time parameter.
- `config::RegularizedCoordinateConfig`: Configuration parameters.

# Returns
- `StiefelScheifele`: The Stiefel-Scheifele state vector.
"""
function cart2StiefelScheifele(
    cart_vec::AbstractVector, μ::Number, ϕ::Number, config::RegularizedCoordinateConfig
)
    DU, TU, W, t₀, flag_time = config.DU, config.TU, config.W, config.t₀, config.flag_time

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
    sph2, cph2 = sincos(ϕ / 2)

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
        t_val = t₀ / TU
    elseif flag_time isa LinearTime
        t_val =
            t₀ / TU + dot(SVector{4}(u₁, u₂, u₃, u₄), SVector{4}(du₁, du₂, du₃, du₄)) / ω
    end

    return SVector{10}(α[1], α[2], α[3], α[4], β[1], β[2], β[3], β[4], ω, t_val)
end

"""
    get_stiefelscheifele_time(u::AbstractVector{T}, ϕ::Number, config::RegularizedCoordinateConfig) where {T<:Number}

Computes the time element for Stiefel-Scheifele transformations.

# Arguments
- `u::AbstractVector{T}`: The Stiefel-Scheifele state vector.
- `ϕ::Number`: Fictitious time parameter.
- `config::RegularizedCoordinateConfig`: Configuration parameters.

# Returns
- The time element.
"""
function get_stiefelscheifele_time(
    u::AbstractVector, ϕ::Number, config::RegularizedCoordinateConfig
)
    t₀, TU, flag_time = config.t₀, config.TU, config.flag_time

    if flag_time isa PhysicalTime
        t = u[10]
    elseif flag_time isa LinearTime
        sph, cph = sincos(ϕ)
        αsq = u[1]^2 + u[2]^2 + u[3]^2 + u[4]^2
        βsq = u[5]^2 + u[6]^2 + u[7]^2 + u[8]^2
        αβ = u[1] * u[5] + u[2] * u[6] + u[3] * u[7] + u[4] * u[8]
        t = u[10] + 0.5 * ((αsq - βsq) / 2.0 * sph - αβ * cph) / u[9]
    end

    return t * TU + t₀
end
