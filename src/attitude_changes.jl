export EP2MRP, MRP2EP

"""
    EP2MRP(β::AbstractVector{<:Number})

Converts Euler Parameter rotation description into Modified Rodriguez Parameters.

# Arguments
- `β::AbstractVector{<:Number}`: The Euler Parameter description of a rotation.

# Returns
- `σ::AbstractVector{<:Number}`: The Modified Rodriguez Parameter description of a rotation.
"""
function EP2MRP(β::AbstractVector{<:Number})
    β0, β1, β2, β3 = β

    # Add small epsilon to avoid division by zero when β0 is close to -1
    denom = 1.0 + β0

    σ1 = β1 / denom
    σ2 = β2 / denom
    σ3 = β3 / denom

    σ = SVector{3}(σ1, σ2, σ3)

    return σ
end

"""
    MRP2EP(σ::AbstractVector{<:Number})

Converts Modified Rodriguez Parameters rotation description into Euler Parameter.

# Arguments
- `σ::AbstractVector{<:Number}`: The Modified Rodriguez Parameter description of a rotation.

# Returns
- `β::AbstractVector{<:Number}`: The Euler Parameter description of a rotation.
"""
function MRP2EP(σ::AbstractVector{<:Number})
    # Avoid broadcasting allocation by computing sum directly
    σ_sq = σ[1]^2 + σ[2]^2 + σ[3]^2

    denom = 1.0 + σ_sq

    β0 = (1.0 - σ_sq) / denom
    β1 = (2.0 * σ[1]) / denom
    β2 = (2.0 * σ[2]) / denom
    β3 = (2.0 * σ[3]) / denom

    return SVector{4}(β0, β1, β2, β3)
end
