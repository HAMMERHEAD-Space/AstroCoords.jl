abstract type KSTimeType end
struct KSPhysicalTime <: KSTimeType end
struct KSLinearTime <: KSTimeType end

export get_KS_time, set_ks_configurations, KSPhysicalTime, KSLinearTime

"""
    cart2KS(u, μ; Vpot, DU, TU, t₀, flag_time=KSPhysicalTime())

Converts a Cartesian state vector to a Kustaanheimo-Stiefel (KS) state vector.

# Arguments
- `u::AbstractVector`: Cartesian state `[x, y, z, ẋ, ẏ, ż]`.
- `μ::Number`: Gravitational parameter.

# Keyword Arguments
- `Vpot::Number`: Perturbing potential. Defaults to zero for the two-body problem.
- `DU::Number`: Distance unit.
- `TU::Number`: Time unit.
- `t₀::Number`: Initial physical time.
- `flag_time::KSTimeType`: Time element formulation (`KSPhysicalTime` or `KSLinearTime`).

# Returns
- `SVector{10, RT}`: The 10-element KS state vector.
"""
function cart2KS(
    u_cart::AbstractVector{T},
    μ::V;
    Vpot::VP,
    DU::DT,
    TU::TT,
    t₀::TT2,
    flag_time::KSTimeType=KSPhysicalTime(),
) where {T<:Number,V<:Number,VP<:Number,DT<:Number,TT<:Number,TT2<:Number}
    RT = promote_type(T, V, VP, DT, TT, TT2)

    ##################################################
    #* 1. Non-Dimensionalize
    ##################################################
    r_vec = SVector{3,RT}(u_cart[1] / DU, u_cart[2] / DU, u_cart[3] / DU)
    v_vec = SVector{3,RT}(u_cart[4] / (DU / TU), u_cart[5] / (DU / TU), u_cart[6] / (DU / TU))
    r = norm(r_vec)
    t = t₀ / TU

    ##################################################
    #* 2. Initialize KS-Position and KS-Velocity
    ##################################################
    # KS-Position
    if r_vec[1] >= 0
        u₁ = 0.0
        u₄ = √(.5 * (r + r_vec[1]))
        u₂ = (r_vec[2]*u₁ + r_vec[3]*u₄) / (r + r_vec[1])
        u₃ = (r_vec[3]*u₁ - r_vec[2]*u₄) / (r + r_vec[1])
    else
        u₂ = 0.0
        u₃ = √(.5 * (r - r_vec[1]))
        u₁ = (r_vec[2]*u₂ + r_vec[3]*u₃) / (r - r_vec[1])
        u₄ = (r_vec[3]*u₂ - r_vec[2]*u₃) / (r - r_vec[1])
    end

    # KS-Velocity
    u₅ = .5 * dot(SVector{3, RT}(u₁, u₂, u₃), v_vec)
    u₆ = .5 * dot(SVector{3, RT}(-u₂, u₁, u₄), v_vec)
    u₇ = .5 * dot(SVector{3, RT}(-u₃, -u₄, u₁), v_vec)
    u₈ = .5 * dot(SVector{3, RT}(u₄, -u₃, u₂), v_vec)

    # Total energy
    h = μ / r - 0.5 * dot(v_vec, v_vec) - Vpot

    # Time element
    if flag_time isa KSPhysicalTime
        τ = t
    elseif flag_time isa KSLinearTime
        τ = t + dot(SVector{4}(u₁, u₂, u₃, u₄), SVector{4}(u₅, u₆, u₇, u₈)) / h
    end

    return SVector{10,RT}(
        u₁, u₂, u₃, u₄,
        u₅, u₆, u₇, u₈,
        h, τ
    )
end

"""
    KS2cart(u_ks, μ)

Converts a Kustaanheimo-Stiefel (KS) state vector to a Cartesian state vector.

# Arguments
- `u_ks::AbstractVector`: 10-element KS state vector.
- `μ::Number`: Gravitational parameter (unused in this direction, kept for API consistency).

# Keyword Arguments
- `Vpot::Number`: Perturbing potential. Defaults to zero for the two-body problem.
- `DU::Number`: Distance unit.
- `TU::Number`: Time unit.
- `t₀::Number`: Initial physical time.
- `flag_time::KSTimeType`: Time element formulation (`KSPhysicalTime` or `KSLinearTime`).

# Returns
- `SVector{6, RT}`: The 6-element Cartesian state `[x, y, z, ẋ, ẏ, ż]`.
"""
function KS2cart(
    u_ks::AbstractVector{T},
    μ::V;
    Vpot::VP,
    DU::DT,
    TU::TT,
    t₀::TT2,
    flag_time::KSTimeType=KSPhysicalTime(),
) where {T<:Number,V<:Number,VP<:Number,DT<:Number,TT<:Number,TT2<:Number}
    RT = promote_type(T, V, VP, DT, TT, TT2)

    # Position
    x = u_ks[1]^2 - u_ks[2]^2 - u_ks[3]^2 + u_ks[4]^2
    y = 2 * (u_ks[1] * u_ks[2] - u_ks[3] * u_ks[4])
    z = 2 * (u_ks[1] * u_ks[3] + u_ks[2] * u_ks[4])
    r = u_ks[1]^2 + u_ks[2]^2 + u_ks[3]^2 + u_ks[4]^2

    # Velocity
    ẋ = 2 * (u_ks[1] * u_ks[5] - u_ks[2] * u_ks[6] - u_ks[3] * u_ks[7] + u_ks[4] * u_ks[8]) / r
    ẏ = 2 * (u_ks[2] * u_ks[5] + u_ks[1] * u_ks[6] - u_ks[4] * u_ks[7] - u_ks[3] * u_ks[8]) / r
    ż = 2 * (u_ks[3] * u_ks[5] + u_ks[4] * u_ks[6] + u_ks[1] * u_ks[7] + u_ks[2] * u_ks[8]) / r
    
    return SVector{6,RT}(x * DU, y * DU, z * DU, ẋ * (DU / TU), ẏ * (DU / TU), ż * (DU / TU))
end

"""
    get_KS_time(u, flag_time)

Computes the physical time from the KS state vector.

# Arguments
- `u::AbstractVector`: KS state vector `[u₁, u₂, u₃, u₄, u₅, u₆, u₇, u₈, h, τ]`.
- `flag_time::KSTimeType`: Time element formulation (`KSPhysicalTime` or `KSLinearTime`).

# Returns
- `Number`: The computed physical time.
"""
function get_KS_time(u::AbstractVector{T}, flag_time::KSTimeType) where {T<:Number}
    if flag_time isa KSPhysicalTime
        t = u[10]
    elseif flag_time isa KSLinearTime
        t = u[10] - dot(SVector{4}(u[1], u[2], u[3], u[4]), SVector{4}(u[5], u[6], u[7], u[8])) / u[9]
    end
    return t
end

"""
    set_ks_configurations(u, μ; Vpot=0.0, t₀=0.0, flag_time=KSPhysicalTime())

Computes and returns a `NamedTuple` of configurations required for KS transformations.

The configurations are derived from the initial Cartesian state vector and the
gravitational parameter.

- `DU` is set to the initial position magnitude.
- `TU` is derived from `DU` and `μ`.


# Arguments
- `u::AbstractVector`: Cartesian state vector `[x, y, z, ẋ, ẏ, ż]`.
- `μ::Number`: Gravitational parameter.

# Keyword Arguments
- `Vpot::Number=0.0`: Perturbing potential.
- `t₀::Number=0.0`: Initial physical time.
- `flag_time::KSTimeType=KSPhysicalTime()`: Time element formulation.

# Returns
- `NamedTuple`: A tuple containing `DU`, `TU`, `Vpot`, `t₀`, `flag_time`.
"""
function set_ks_configurations(u::AbstractVector, μ; Vpot=0.0, t₀=0.0, flag_time=KSPhysicalTime())
    DU = norm(SVector{3}(u[1], u[2], u[3]))
    TU = sqrt(DU^3 / μ)
    return (; DU, TU, Vpot, t₀, flag_time)
end

