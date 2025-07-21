export RegularizedCoordinateConfig
export compute_characteristic_scales, compute_initial_phi

"""
    RegularizedCoordinateConfig{DT,TT,WT,PT,TT2,FT}

Configuration struct for regularized coordinate transformations (EDromo, Kustaanheimo-Stiefel, Stiefel-Scheifele).

This struct encapsulates all the parameters required for regularized coordinate transformations:
- `DU`: Reference distance unit for non-dimensionalization
- `TU`: Reference time unit for non-dimensionalization  
- `W`: Perturbing potential energy (unified for all coordinate systems)
- `ϕ`: Initial value of the fictitious time (used by EDromo and SS, ignored by KS)
- `t₀`: Initial physical time
- `flag_time`: Time element formulation (`PhysicalTime`, `ConstantTime`, or `LinearTime`)
"""
struct RegularizedCoordinateConfig{
    DT<:Number,TT<:Number,WT<:Number,PT<:Number,TT2<:Number,FT<:AbstractTimeType
}
    DU::DT
    TU::TT
    W::WT
    ϕ::PT
    t₀::TT2
    flag_time::FT
end

# Keyword constructor
function RegularizedCoordinateConfig(;
    DU::Number=0.0,
    TU::Number=0.0,
    W::Number=0.0,
    ϕ::Number=0.0,
    t₀::Number=0.0,
    flag_time::AbstractTimeType=PhysicalTime(),
)
    return RegularizedCoordinateConfig(DU, TU, W, ϕ, t₀, flag_time)
end

# Method to update ϕ while keeping other parameters
function RegularizedCoordinateConfig(config::RegularizedCoordinateConfig; curr_ϕ::Number)
    return RegularizedCoordinateConfig(
        config.DU, config.TU, config.W, curr_ϕ, config.t₀, config.flag_time
    )
end

"""
    compute_characteristic_scales(state, μ)

Computes characteristic distance and time scales for regularized coordinate transformations.

The scales are computed as:
- `DU`: Distance unit based on the position magnitude  
- `TU`: Time unit derived from `DU` and `μ` using Kepler's third law: `TU = sqrt(DU³/μ)`

# Arguments
- `state::AbstractVector`: Cartesian state vector `[x, y, z, ẋ, ẏ, ż]`
- `μ::Number`: Gravitational parameter

# Returns
- `(DU, TU)`: Tuple of characteristic distance and time scales
"""
function compute_characteristic_scales(state::AbstractVector, μ::Number)
    DU = norm(SVector{3}(state[1], state[2], state[3]))
    TU = sqrt(DU^3 / μ)
    return (DU, TU)
end

"""
    compute_initial_phi(state, μ, DU, TU, W)

Computes the initial value of the fictitious time `ϕ` for regularized coordinates (EDromo, Stiefel-Scheifele).

This calculation is based on Baù, G., Bombardelli, C., Peláez, J., and Lorenzini, E., 
"Nonsingular orbital elements for special perturbations in the two-body problem". 
MNRAS 454(3), pp. 2890-2908. 2015.

# Arguments
- `state::AbstractVector`: Cartesian state vector `[x, y, z, ẋ, ẏ, ż]`
- `μ::Number`: Gravitational parameter
- `DU::Number`: Reference distance unit
- `TU::Number`: Reference time unit  
- `W::Number`: Perturbing potential energy

# Returns
- `Number`: The computed value of `ϕ₀`
"""
function compute_initial_phi(
    state::AbstractVector{T}, μ::V, DU::DT, TU::TT, W::WT
) where {T<:Number,V<:Number,DT<:Number,TT<:Number,WT<:Number}
    RT = promote_type(T, V, DT, TT, WT)

    x, y, z, ẋ, ẏ, ż = state

    r = SVector{3,RT}(x / DU, y / DU, z / DU)
    v = SVector{3,RT}(ẋ / (DU / TU), ẏ / (DU / TU), ż / (DU / TU))

    Μ = μ / (DU^3 / TU^2)
    W_nd = W / (DU^2 / TU^2)

    E = 0.5 * dot(v, v) - Μ / norm(r) - W_nd

    return atan(dot(r, v)*√(-2.0*E), 1.0 + 2.0*E*norm(r))
end

"""
    RegularizedCoordinateConfig(state, μ; W=0.0, ϕ=0.0, t₀=0.0, flag_time=PhysicalTime())

Convenience constructor that computes characteristic scales from initial state and creates configuration.

# Arguments
- `state::AbstractVector`: Cartesian state vector `[x, y, z, ẋ, ẏ, ż]`
- `μ::Number`: Gravitational parameter

# Keyword Arguments  
- `W::Number=0.0`: Perturbing potential energy
- `ϕ::Number=0.0`: Initial fictitious time (use `compute_initial_phi` for EDromo/SS)
- `t₀::Number=0.0`: Initial physical time
- `flag_time::AbstractTimeType=PhysicalTime()`: Time element formulation

# Returns
- `RegularizedCoordinateConfig`: Configuration struct with computed scales
"""
function RegularizedCoordinateConfig(
    state::AbstractVector,
    μ::Number;
    W::Number=0.0,
    ϕ::Union{Number,Nothing}=nothing,
    t₀::Number=0.0,
    flag_time::AbstractTimeType=PhysicalTime(),
)
    DU, TU = compute_characteristic_scales(state, μ)

    # Auto-compute ϕ if not provided
    if ϕ === nothing
        ϕ = compute_initial_phi(state, μ, DU, TU, W)
    end

    return RegularizedCoordinateConfig(; DU=DU, TU=TU, W=W, ϕ=ϕ, t₀=t₀, flag_time=flag_time)
end
