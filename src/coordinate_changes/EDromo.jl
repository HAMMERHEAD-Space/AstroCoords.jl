export get_EDromo_time

"""
    cart2EDromo(u, μ, config::RegularizedCoordinateConfig)

Converts a Cartesian state vector to an EDromo state vector.

This is the backend implementation for the transformation from `Cartesian` to `EDromo`.
It requires a `RegularizedCoordinateConfig` with all necessary parameters.

# Arguments
- `u::AbstractVector`: Cartesian state vector `[x, y, z, ẋ, ẏ, ż]` in `[L]` and `[L/T]`.
- `μ::Number`: Gravitational parameter in `[L³/T²]`.
- `config::RegularizedCoordinateConfig`: Configuration parameters for the transformation.

# Returns
- `SVector{8, RT}`: The 8-element EDromo state vector.
"""
function cart2EDromo(
    u::AbstractVector{T}, μ::V, config::RegularizedCoordinateConfig
) where {T<:Number,V<:Number}
    DU, TU, W, ϕ, t₀, flag_time = config.DU,
    config.TU, config.W, config.ϕ, config.t₀,
    config.flag_time
    RT = promote_type(T, V, typeof(DU), typeof(TU), typeof(W), typeof(ϕ), typeof(t₀))

    x, y, z, ẋ, ẏ, ż = u

    ##################################################
    #* 1. Non-Dimensionalize
    ##################################################

    r = SVector{3,RT}(x / DU, y / DU, z / DU)
    v = SVector{3,RT}(ẋ / (DU / TU), ẏ / (DU / TU), ż / (DU / TU))

    W_non_dim = W / (DU^2 / TU^2)
    Μ = μ / (DU^3 / TU^2)

    ##################################################
    #* 2. In-plane Elements
    ##################################################

    r_mag = norm(r)
    v_mag = norm(v)

    sϕ, cϕ = sincos(ϕ)

    # Total Energy
    E = v_mag^2 / 2 - Μ / r_mag + W_non_dim

    # Angular Momentum
    h₀ = cross(r, v)
    h₀_mag = norm(h₀)

    # Generalized Angular Momentum
    c₀ = √(h₀_mag^2 + 2*r_mag^2*W_non_dim)

    # Dot Product of r and v
    r_dot_v = dot(r, v)

    ζ₁ = (1.0 + 2.0*E*r_mag)*cϕ + r_dot_v * √(-2.0 * E) * sϕ
    ζ₂ = (1.0 + 2.0*E*r_mag)*sϕ - r_dot_v * √(-2.0 * E) * cϕ
    ζ₃ = -1.0 / (2.0 * E)

    ##################################################
    #* 3. Quaternion Elements
    ##################################################
    ν₀ = ϕ + 2.0*atan(r_dot_v, c₀ + r_mag*√(-2.0*E))
    sν, cν = sincos(ν₀)

    # Orbital Frame in the IRF
    î = r / r_mag
    k̂ = h₀ / h₀_mag
    ĵ = cross(k̂, î)

    # Intermediate Frame Unit Vectors in the IRF
    x̂ = cν * î - sν * ĵ
    ŷ = sν * î + cν * ĵ

    # Safe Initialization of the Quaternion
    # Arguments of the roots have to be always positive, so we can safely use abs()
    aux = abs(1.0 + x̂[1] + ŷ[2] + k̂[3])
    ζ₇ = 0.5 * √(aux)

    # Check for Singularities and NaNs
    if aux <= eps(RT)
        aux = abs(0.5 * (k̂[3] + 1.0))
        ζ₆ = √(aux)
        if aux <= eps(RT)
            aux = abs(0.5 * (1.0 - ŷ[2]))
            ζ₄ = √(aux)
            if aux <= eps(RT)
                ζ₅ = 1.0
            else
                ζ₅ = ŷ[1] / (2.0*ζ₄)
            end
        else
            ζ₄ = k̂[1] / (2.0*ζ₆)
            ζ₅ = k̂[2] / (2.0*ζ₆)
        end
    else
        ζ₄ = (ŷ[3] - k̂[2]) / (4.0 * ζ₇)
        ζ₅ = (k̂[1] - x̂[3]) / (4.0 * ζ₇)
        ζ₆ = (x̂[2] - ŷ[1]) / (4.0 * ζ₇)
    end

    ##################################################
    #* 4. Time Element
    ##################################################
    if flag_time isa PhysicalTime
        ζ₈ = t₀ / TU
    elseif flag_time isa ConstantTime
        ζ₈ = t₀ / TU + ζ₃^(1.5) * (ζ₁*sϕ - ζ₂*cϕ - ϕ)
    elseif flag_time isa LinearTime
        ζ₈ = t₀ / TU + ζ₃^(1.5) * (ζ₁*sϕ - ζ₂*cϕ)
    end

    return SVector{8,RT}(ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈)
end

"""
    EDromo2cart(u, μ, config::RegularizedCoordinateConfig)

Converts an EDromo state vector to a Cartesian state vector.

This is the backend implementation for the transformation from `EDromo` to `Cartesian`.
It requires a `RegularizedCoordinateConfig` with the same parameters used in the forward transformation.

# Arguments
- `u::AbstractVector`: EDromo state vector `[ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈]`.
- `μ::Number`: Gravitational parameter of the central body.
- `config::RegularizedCoordinateConfig`: Configuration parameters for the transformation.

# Returns
- `SVector{6, RT}`: The 6-element Cartesian state vector `[x, y, z, ẋ, ẏ, ż]`.
"""
function EDromo2cart(
    u::AbstractVector{T}, μ::V, config::RegularizedCoordinateConfig
) where {T<:Number,V<:Number}
    DU, TU, W, ϕ, t₀ = config.DU, config.TU, config.W, config.ϕ, config.t₀
    RT = promote_type(T, V, typeof(DU), typeof(TU), typeof(W), typeof(ϕ), typeof(t₀))

    ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈ = u

    ##################################################
    #* 1. Auxiliary Quantities
    ##################################################
    sϕ, cϕ = sincos(ϕ)

    ρ = 1.0 - ζ₁*cϕ - ζ₂*sϕ
    r_mag = ζ₃ * ρ

    Z = ζ₁*sϕ - ζ₂*cϕ
    emme = √(1.0 - ζ₁^2 - ζ₂^2)

    cν = (cϕ - ζ₁ + (Z*ζ₂) / (emme + 1.0)) / ρ
    sν = (sϕ - ζ₂ - (Z*ζ₁) / (emme + 1.0)) / ρ

    ##################################################
    #* 2. Compute Position in the Inertial Frame
    ##################################################

    # Intermediate Frame Unit Vectors
    x̂ = 2.0 * SVector{3}(0.5 - ζ₅^2 - ζ₆^2, ζ₄*ζ₅ + ζ₆*ζ₇, ζ₄*ζ₆ - ζ₅*ζ₇)

    ŷ = 2.0 * SVector{3}(ζ₄*ζ₅ - ζ₆*ζ₇, 0.5 - ζ₄^2 - ζ₆^2, ζ₅*ζ₆ + ζ₄*ζ₇)

    r_non_dim = r_mag * (x̂ * cν + ŷ * sν)

    ##################################################
    #* 3. Perturbing Potential
    ##################################################
    U = W / (DU^2 / TU^2)

    ##################################################
    #* 4. Compute Velocity in the Inertial Frame
    ##################################################
    î = x̂*cν + ŷ*sν
    ĵ = -x̂*sν + ŷ*cν

    v_rad = Z / (√(ζ₃)*ρ)
    v_tan = √((1.0 - ζ₁^2 - ζ₂^2) / (ζ₃*ρ^2) - 2.0 * U)

    v_non_dim = v_rad * î + v_tan * ĵ

    r = r_non_dim * DU
    v = v_non_dim * (DU / TU)

    return SVector{6,RT}(r[1], r[2], r[3], v[1], v[2], v[3])
end

"""
    get_EDromo_time(u, config::RegularizedCoordinateConfig)

Computes the physical time from the EDromo state vector.

# Arguments
- `u::AbstractVector`: EDromo state vector `[ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈]`.
- `config::RegularizedCoordinateConfig`: Configuration parameters for the transformation.

# Returns
- `Number`: The computed physical time.
"""
function get_EDromo_time(
    u::AbstractVector{T}, config::RegularizedCoordinateConfig
) where {T<:Number}
    TU, ϕ, t₀, flag_time = config.TU, config.ϕ, config.t₀, config.flag_time
    RT = promote_type(T, typeof(ϕ))

    sϕ, cϕ = sincos(ϕ)

    if flag_time isa PhysicalTime
        t = u[8]
    elseif flag_time isa ConstantTime
        t = u[8] - u[3]^(1.5) * (u[1]*sϕ - u[2]*cϕ - ϕ)
    elseif flag_time isa LinearTime
        t = u[8] - u[3]^(1.5) * (u[1]*sϕ - u[2]*cϕ)
    end

    return RT(t * TU + t₀)
end
