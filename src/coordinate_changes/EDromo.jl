abstract type AbstractTimeType end
struct PhysicalTime <: AbstractTimeType end
struct ConstantTime <: AbstractTimeType end
struct LinearTime <: AbstractTimeType end

export get_EDromo_time, set_edromo_configurations, PhysicalTime, ConstantTime, LinearTime

"""
    cart2EDromo(u, μ; DU, TU, ϕ₀, t₀, W, flag_time)

Converts a Cartesian state vector to an EDromo state vector.

This is the backend implementation for the transformation from `Cartesian` to `EDromo`.
It requires a set of non-dimensionalizing parameters and configuration flags.

# Arguments
- `u::AbstractVector`: Cartesian state vector `[x, y, z, ẋ, ẏ, ż]` in `[L]` and `[L/T]`.
- `μ::Number`: Gravitational parameter in `[L³/T²]`.

# Keyword Arguments
- `DU::Number`: Reference distance unit for non-dimensionalization.
- `TU::Number`: Reference time unit for non-dimensionalization.
- `ϕ₀::Number`: Initial value of the fictitious time `ϕ`.
- `t₀::Number`: Initial physical time `t`.
- `W::Number`: Perturbing potential energy.
- `flag_time::AbstractTimeType`: Time element formulation (`PhysicalTime`, `ConstantTime`, or `LinearTime`).

# Returns
- `SVector{8, RT}`: The 8-element EDromo state vector.
"""
function cart2EDromo(
    u::AbstractVector{T},
    μ::V;
    DU::DT,
    TU::TT,
    ϕ₀::PT,
    t₀::TT2,
    W::WT,
    flag_time::AbstractTimeType,
) where {T<:Number,V<:Number,DT<:Number,TT<:Number,PT<:Number,TT2<:Number,WT<:Number}
    RT = promote_type(T, V, DT, TT, PT, TT2, WT)

    x, y, z, ẋ, ẏ, ż = u

    ##################################################
    #* 1. Non-Dimensionalize
    ##################################################

    r = SVector{3,RT}(x / DU, y / DU, z / DU)
    v = SVector{3,RT}(ẋ / (DU / TU), ẏ / (DU / TU), ż / (DU / TU))

    W₀ = W / (DU / TU)^2

    ##################################################
    #* 2. In-plane Elements
    ##################################################

    r_mag = norm(r)
    v_mag = norm(v)

    sϕ, cϕ = sincos(ϕ₀)

    # Total Energy
    E = v_mag^2 / 2 - 1.0 / r_mag + W₀

    # Angular Momentum
    h₀ = cross(r, v)
    h₀_mag = norm(h₀)

    # Generalized Angular Momentum
    c₀ = √(h₀_mag^2 + 2*r_mag^2*W₀)

    # Dot Product of r and v
    r_dot_v = dot(r, v)

    ζ₁ = (1.0 + 2.0*E*r_mag)*cϕ + r_dot_v * √(-2.0 * E) * sϕ
    ζ₂ = (1.0 + 2.0*E*r_mag)*sϕ - r_dot_v * √(-2.0 * E) * cϕ
    ζ₃ = -1.0 / (2.0 * E)

    ##################################################
    #* 3. Quaternion Elements
    ##################################################
    ν₀ = ϕ₀ + 2.0*atan(r_dot_v, c₀ + r_mag*√(-2.0*E))
    sν, cν = sincos(ν₀)

    # Orbital Frame in the IRF
    î = r / r_mag
    k̂ = h₀ / h₀_mag
    ĵ = cross(k̂, î)

    # Intermediate Frame Unit Vectors in the IRF
    x̂ = cν * î - sν * ĵ
    ŷ = sν * î + cν * ĵ

    # Safe Initialization of the Quaternion
    # Arguments of the roots have to be always positive, so we can safely use abs()
    aux = abs(1.0 + x̂[1] + ŷ[2] + k̂[3])
    ζ₇ = 0.5 * √(aux)

    # Check for Singularities and NaNs
    if aux <= eps(RT)
        aux = abs(0.5 * (k̂[3] + 1.0))
        ζ₆ = √(aux)
        if aux <= eps(RT)
            aux = abs(0.5 * (1.0 - ŷ[2]))
            ζ₄ = √(aux)
            if aux <= eps(RT)
                ζ₅ = 1.0
            else
                ζ₅ = ŷ[1] / (2.0*ζ₄)
            end
        else
            ζ₄ = k̂[1] / (2.0*ζ₆)
            ζ₅ = k̂[2] / (2.0*ζ₆)
        end
    else
        ζ₄ = (ŷ[3] - k̂[2]) / (4.0 * ζ₇)
        ζ₅ = (k̂[1] - x̂[3]) / (4.0 * ζ₇)
        ζ₆ = (x̂[2] - ŷ[1]) / (4.0 * ζ₇)
    end

    ##################################################
    #* 4. Time Element
    ##################################################
    if flag_time isa PhysicalTime
        ζ₈ = t₀ / TU
    elseif flag_time isa ConstantTime
        ζ₈ = t₀ / TU + ζ₃^(1.5) * (ζ₁*sϕ - ζ₂*cϕ - ϕ₀)
    elseif flag_time isa LinearTime
        ζ₈ = t₀ / TU + ζ₃^(1.5) * (ζ₁*sϕ - ζ₂*cϕ)
    end

    return SVector{8,RT}(ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈)
end

"""
    EDromo2cart(u, μ; DU, TU, ϕ₀, t₀, W, flag_time)

Converts an EDromo state vector to a Cartesian state vector.

This is the backend implementation for the transformation from `EDromo` to `Cartesian`.
It requires the same set of non-dimensionalizing parameters and configuration flags
used in the forward transformation.

# Arguments
- `u::AbstractVector`: EDromo state vector `[ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈]`.
- `μ::Number`: Gravitational parameter of the central body.

# Keyword Arguments
- `DU::Number`: Reference distance unit for non-dimensionalization.
- `TU::Number`: Reference time unit for non-dimensionalization.
- `ϕ₀::Number`: Initial value of the fictitious time `ϕ`.
- `t₀::Number`: Initial physical time `t`.
- `W::Number`: Perturbing potential energy.
- `flag_time::AbstractTimeType`: Time element formulation (`PhysicalTime`, `ConstantTime`, or `LinearTime`).

# Returns
- `SVector{6, RT}`: The 6-element Cartesian state vector `[x, y, z, ẋ, ẏ, ż]`.
"""
function EDromo2cart(
    u::AbstractVector{T},
    μ::V;
    DU::DT,
    TU::TT,
    ϕ₀::PT,
    t₀::TT2,
    W::WT,
    flag_time::AbstractTimeType,
) where {T<:Number,V<:Number,DT<:Number,TT<:Number,PT<:Number,TT2<:Number,WT<:Number}
    RT = promote_type(T, V, DT, TT, PT, TT2, WT)

    ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈ = u

    ##################################################
    #* 1. Auxiliary Quantities
    ##################################################
    sϕ, cϕ = sincos(ϕ₀)

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

    ŷ = 2.0 * SVector{3}(ζ₄*ζ₅ - ζ₆*ζ₇, 0.5 - ζ₄^2 - ζ₆^2, ζ₅*ζ₆ + ζ₄*ζ₇)

    r_non_dim = r_mag * (x̂ * cν + ŷ * sν)

    ##################################################
    #* 3. Perturbing Potential
    ##################################################
    U = W / (DU / TU)^2

    ##################################################
    #* 4. Compute Velocity in the Intertial Frame
    ##################################################
    î = x̂*cν + ŷ*sν
    ĵ = -x̂*sν + ŷ*cν

    v_rad = Z / (√(ζ₃)*ρ)
    v_tan = √((1.0 - ζ₁^2 - ζ₂^2) / (ζ₃*ρ^2) - 2.0 * U)

    v_non_dim = v_rad * î + v_tan * ĵ

    r = r_non_dim * DU
    v = v_non_dim * (DU / TU)

    return SVector{6,RT}(r[1], r[2], r[3], v[1], v[2], v[3])
end

"""
    EDromo_get_time(u, flag_time)

Computes the physical time from the EDromo state vector.

# Arguments
- `u::AbstractVector`: EDromo state vector `[ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈]`.
- `flag_time::AbstractTimeType`: Time element formulation (`PhysicalTime`, `ConstantTime`, or `LinearTime`).

# Returns
- `Number`: The computed physical time.
"""
function get_EDromo_time(
    u::AbstractVector{T}, flag_time::AbstractTimeType
) where {T<:Number}
    if flag_time isa PhysicalTime
        t = u[8]
    elseif flag_time isa ConstantTime
        t = u[8] - u[3]^(1.5) * (u[1]*sϕ - u[2]*cϕ - ϕ₀)
    elseif flag_time isa LinearTime
        t = u[8] - u[3]^(1.5) * (u[1]*sϕ - u[2]*cϕ)
    end

    return t
end

"""
    set_edromo_configurations(u, μ; W=0.0, t₀=0.0, flag_time=PhysicalTime())

Computes and returns a `NamedTuple` of configurations required for EDromo transformations.

The configurations are derived from the initial Cartesian state vector and the
gravitational parameter.

- `DU` is set to the initial position magnitude.
- `TU` is derived from `DU` and `μ`.
- `ϕ₀` is computed based on the formulation in Baù, et al. (2015).

# Arguments
- `u::AbstractVector`: Cartesian state vector `[x, y, z, ẋ, ẏ, ż]`.
- `μ::Number`: Gravitational parameter.

# Keyword Arguments
- `W::Number=0.0`: Perturbing potential energy.
- `t₀::Number=0.0`: Initial physical time.
- `flag_time::AbstractTimeType=PhysicalTime()`: Time element formulation.

# Returns
- `NamedTuple`: A tuple containing `DU`, `TU`, `W`, `ϕ₀`, `t₀`, `flag_time`.
"""
function set_edromo_configurations(
    u::AbstractVector, μ; W=0.0, t₀=0.0, flag_time=PhysicalTime()
)
    DU = norm(SVector{3}(u[1], u[2], u[3]))
    TU = sqrt(DU^3 / μ)
    ϕ₀ = computeϕ₀(u, μ, DU, TU, W)
    return (; DU, TU, W, ϕ₀, t₀, flag_time)
end

"""
    computeϕ₀(u, μ, DU, TU, W₀)

Computes the initial value of the fictitious time `ϕ₀` for the EDromo formulation.

This calculation is based on the suggestion in [1] from the Fortran source, which
corresponds to Baù, G., Bombardelli, C., Peláez, J., and Lorenzini, E., "Nonsingular
orbital elements for special perturbations in the two-body problem". MNRAS 454(3),
pp. 2890-2908. 2015.

# Arguments
- `u::AbstractVector`: Cartesian state vector.
- `μ::Number`: Gravitational parameter.
- `DU::Number`: Reference distance unit.
- `TU::Number`: Reference time unit.
- `W₀::Number`: Perturbing potential energy.

# Returns
- `Number`: The computed value of `ϕ₀`.
"""
function computeϕ₀(
    u::AbstractVector{T}, μ::V, DU::DT, TU::TT, W₀::WT
) where {T<:Number,V<:Number,DT<:Number,TT<:Number,WT<:Number}
    RT = promote_type(T, V, DT, TT, WT)

    x, y, z, ẋ, ẏ, ż = u

    r = SVector{3,RT}(x / DU, y / DU, z / DU)
    v = SVector{3,RT}(ẋ / (DU / TU), ẏ / (DU / TU), ż / (DU / TU))

    Μ = μ / (DU^3 / TU^2)
    U = W₀ / (DU^2 / TU^2)

    E = 0.5 * dot(v, v) - Μ / norm(r) - U

    return atan(dot(r, v)*√(-2.0*E), 1.0 + 2.0*E*norm(r))
end
