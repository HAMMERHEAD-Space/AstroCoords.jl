export get_KS_time

"""
    cart2KS(u, μ, config::RegularizedCoordinateConfig)

Converts a Cartesian state vector to a Kustaanheimo-Stiefel (KS) state vector.

# Arguments
- `u::AbstractVector`: Cartesian state `[x, y, z, ẋ, ẏ, ż]`.
- `μ::Number`: Gravitational parameter.
- `config::RegularizedCoordinateConfig`: Configuration parameters (uses config.W for perturbing potential).

# Returns
- `SVector{10, RT}`: The 10-element KS state vector.
"""
function cart2KS(
    u_cart::AbstractVector{T}, μ::V, config::RegularizedCoordinateConfig
) where {T<:Number,V<:Number}
    W, DU, TU, t₀, flag_time = config.W, config.DU, config.TU, config.t₀, config.flag_time
    RT = promote_type(T, V, typeof(W), typeof(DU), typeof(TU), typeof(t₀))

    ##################################################
    #* 1. Non-Dimensionalize
    ##################################################
    r_vec = SVector{3,RT}(u_cart[1] / DU, u_cart[2] / DU, u_cart[3] / DU)
    v_vec = SVector{3,RT}(
        u_cart[4] / (DU / TU), u_cart[5] / (DU / TU), u_cart[6] / (DU / TU)
    )
    r = norm(r_vec)
    t = t₀ / TU

    ##################################################
    #* 2. Initialize KS-Position and KS-Velocity
    ##################################################
    # KS-Position
    if r_vec[1] >= 0
        u₁ = 0.0
        u₄ = √(0.5 * (r + r_vec[1]))
        u₂ = (r_vec[2]*u₁ + r_vec[3]*u₄) / (r + r_vec[1])
        u₃ = (r_vec[3]*u₁ - r_vec[2]*u₄) / (r + r_vec[1])
    else
        u₂ = 0.0
        u₃ = √(0.5 * (r - r_vec[1]))
        u₁ = (r_vec[2]*u₂ + r_vec[3]*u₃) / (r - r_vec[1])
        u₄ = (r_vec[3]*u₂ - r_vec[2]*u₃) / (r - r_vec[1])
    end

    # KS-Velocity
    u₅ = 0.5 * dot(SVector{3,RT}(u₁, u₂, u₃), v_vec)
    u₆ = 0.5 * dot(SVector{3,RT}(-u₂, u₁, u₄), v_vec)
    u₇ = 0.5 * dot(SVector{3,RT}(-u₃, -u₄, u₁), v_vec)
    u₈ = 0.5 * dot(SVector{3,RT}(u₄, -u₃, u₂), v_vec)

    # Total energy
    h = μ / r - 0.5 * dot(v_vec, v_vec) - W

    # Time element
    if flag_time isa PhysicalTime
        τ = t
    elseif flag_time isa LinearTime
        τ = t + dot(SVector{4}(u₁, u₂, u₃, u₄), SVector{4}(u₅, u₆, u₇, u₈)) / h
    end

    return SVector{10,RT}(u₁, u₂, u₃, u₄, u₅, u₆, u₇, u₈, h, τ)
end

"""
    KS2cart(u_ks, μ, config::RegularizedCoordinateConfig)

Converts a Kustaanheimo-Stiefel (KS) state vector to a Cartesian state vector.

# Arguments
- `u_ks::AbstractVector`: 10-element KS state vector.
- `μ::Number`: Gravitational parameter (unused in this direction, kept for API consistency).
- `config::RegularizedCoordinateConfig`: Configuration parameters.

# Returns
- `SVector{6, RT}`: The 6-element Cartesian state `[x, y, z, ẋ, ẏ, ż]`.
"""
function KS2cart(
    u_ks::AbstractVector{T}, μ::V, config::RegularizedCoordinateConfig
) where {T<:Number,V<:Number}
    W, DU, TU, t₀, flag_time = config.W, config.DU, config.TU, config.t₀, config.flag_time
    RT = promote_type(T, V, typeof(W), typeof(DU), typeof(TU), typeof(t₀))

    # Position
    x = u_ks[1]^2 - u_ks[2]^2 - u_ks[3]^2 + u_ks[4]^2
    y = 2 * (u_ks[1] * u_ks[2] - u_ks[3] * u_ks[4])
    z = 2 * (u_ks[1] * u_ks[3] + u_ks[2] * u_ks[4])
    r = u_ks[1]^2 + u_ks[2]^2 + u_ks[3]^2 + u_ks[4]^2

    # Velocity
    ẋ =
        2 *
        (u_ks[1] * u_ks[5] - u_ks[2] * u_ks[6] - u_ks[3] * u_ks[7] + u_ks[4] * u_ks[8]) / r
    ẏ =
        2 *
        (u_ks[2] * u_ks[5] + u_ks[1] * u_ks[6] - u_ks[4] * u_ks[7] - u_ks[3] * u_ks[8]) / r
    ż =
        2 *
        (u_ks[3] * u_ks[5] + u_ks[4] * u_ks[6] + u_ks[1] * u_ks[7] + u_ks[2] * u_ks[8]) / r

    return SVector{6,RT}(
        x * DU, y * DU, z * DU, ẋ * (DU / TU), ẏ * (DU / TU), ż * (DU / TU)
    )
end

"""
    get_KS_time(u, config::RegularizedCoordinateConfig)

Computes the physical time from the KS state vector.

# Arguments
- `u::AbstractVector`: KS state vector `[u₁, u₂, u₃, u₄, u₅, u₆, u₇, u₈, h, τ]`.
- `config::RegularizedCoordinateConfig`: Configuration parameters.

# Returns
- `Number`: The computed physical time.
"""
function get_KS_time(
    u::AbstractVector{T}, config::RegularizedCoordinateConfig
) where {T<:Number}
    flag_time = config.flag_time
    if flag_time isa PhysicalTime
        t = u[10]
    elseif flag_time isa LinearTime
        t =
            u[10] -
            dot(SVector{4}(u[1], u[2], u[3], u[4]), SVector{4}(u[5], u[6], u[7], u[8])) /
            u[9]
    end
    return t
end
