module AstroCoordsZygoteExt

using AstroCoords
using LinearAlgebra

using Zygote.ChainRulesCore: ChainRulesCore
import Zygote.ChainRulesCore: Tangent, NoTangent, ProjectTo
import Zygote: jacobian

function ChainRulesCore.rrule(
    new_coord::Type{<:AstroCoords.AstroCoord}, coord::AbstractArray
)
    # Forward evaluation (basic coordinate transformation)
    y = new_coord(coord)

    function AstroCoord_pullback(Δ::AbstractVector)
        # Define the pullback (how gradients propagate backwards)
        Δcoords = typeof(coord)(Δ)
        return (NoTangent(), Δcoords)
    end

    function AstroCoord_pullback(Δ::Tangent)
        return (NoTangent(), collect(values(Δ)))
    end

    return y, AstroCoord_pullback
end

# rrule for EDromo transformations
function ChainRulesCore.rrule(
    ::Type{EDromo}, coord::Cartesian, μ::Real; kwargs...
)
    # Forward evaluation
    y = EDromo(coord, μ; kwargs...)
    
    function EDromo_pullback(Δy)
        # Extract parameters for the transformation
        coord_params = AstroCoords.params(coord)
        
        # Convert Tangent to vector if needed
        if Δy isa Tangent
            Δy_vec = collect(values(Δy))
        else
            Δy_vec = Δy
        end
        
        # Compute the Jacobian ∂(EDromo)/∂(Cartesian) analytically
        # Use Zygote to differentiate the underlying transformation function
        J = jacobian(x -> AstroCoords.cart2EDromo(x, μ; kwargs...), coord_params)
        
        # Apply chain rule: Δcoord = J^T * Δy
        Δcoord_params = J' * Δy_vec
        
        return (NoTangent(), Δcoord_params, NoTangent())
    end

    return y, EDromo_pullback
end

function ChainRulesCore.rrule(
    ::Type{Cartesian}, coord::EDromo, μ::Real; kwargs...
)
    # Forward evaluation
    y = Cartesian(coord, μ; kwargs...)
    
    function Cartesian_from_EDromo_pullback(Δy)
        # Extract parameters for the transformation
        coord_params = AstroCoords.params(coord)
        
        # Convert Tangent to vector if needed
        if Δy isa Tangent
            Δy_vec = collect(values(Δy))
        else
            Δy_vec = Δy
        end
        
        # Compute the Jacobian ∂(Cartesian)/∂(EDromo) analytically
        # Use Zygote to differentiate the underlying transformation function
        J = jacobian(x -> AstroCoords.EDromo2cart(x, μ; kwargs...), coord_params)
        
        # Apply chain rule: Δcoord = J^T * Δy
        Δcoord_params = J' * Δy_vec
        
        return (NoTangent(), Δcoord_params, NoTangent())
    end

    return y, Cartesian_from_EDromo_pullback
end

# rrule for KustaanheimoStiefel transformations
function ChainRulesCore.rrule(
    ::Type{KustaanheimoStiefel}, coord::Cartesian, μ::Real; kwargs...
)
    # Forward evaluation
    y = KustaanheimoStiefel(coord, μ; kwargs...)
    
    function KS_pullback(Δy)
        # Extract parameters for the transformation
        coord_params = AstroCoords.params(coord)
        
        # Convert Tangent to vector if needed
        if Δy isa Tangent
            Δy_vec = collect(values(Δy))
        else
            Δy_vec = Δy
        end
        
        # Compute the Jacobian ∂(KustaanheimoStiefel)/∂(Cartesian) analytically
        # Use Zygote to differentiate the underlying transformation function
        J = jacobian(x -> AstroCoords.cart2KS(x, μ; kwargs...), coord_params)
        
        # Apply chain rule: Δcoord = J^T * Δy
        Δcoord_params = J' * Δy_vec
        
        return (NoTangent(), Δcoord_params, NoTangent())
    end

    return y, KS_pullback
end

function ChainRulesCore.rrule(
    ::Type{Cartesian}, coord::KustaanheimoStiefel, μ::Real; kwargs...
)
    # Forward evaluation
    y = Cartesian(coord, μ; kwargs...)
    
    function Cartesian_from_KS_pullback(Δy)
        # Extract parameters for the transformation
        coord_params = AstroCoords.params(coord)
        
        # Convert Tangent to vector if needed
        if Δy isa Tangent
            Δy_vec = collect(values(Δy))
        else
            Δy_vec = Δy
        end
        
        # Compute the Jacobian ∂(Cartesian)/∂(KustaanheimoStiefel) analytically
        # Use Zygote to differentiate the underlying transformation function
        J = jacobian(x -> AstroCoords.KS2cart(x, μ; kwargs...), coord_params)
        
        # Apply chain rule: Δcoord = J^T * Δy
        Δcoord_params = J' * Δy_vec
        
        return (NoTangent(), Δcoord_params, NoTangent())
    end

    return y, Cartesian_from_KS_pullback
end

# rrule for StiefelScheifele transformations
function ChainRulesCore.rrule(
    ::Type{StiefelScheifele}, coord::Cartesian, μ::Real; kwargs...
)
    # Forward evaluation
    y = StiefelScheifele(coord, μ; kwargs...)
    
    function SS_pullback(Δy)
        # Extract parameters for the transformation
        coord_params = AstroCoords.params(coord)
        
        # Convert Tangent to vector if needed
        if Δy isa Tangent
            Δy_vec = collect(values(Δy))
        else
            Δy_vec = Δy
        end
        
        # Compute the Jacobian ∂(StiefelScheifele)/∂(Cartesian) analytically
        # Use Zygote to differentiate the underlying transformation function
        J = jacobian(x -> AstroCoords.cart2StiefelScheifele(x, μ; kwargs...), coord_params)
        
        # Apply chain rule: Δcoord = J^T * Δy
        Δcoord_params = J' * Δy_vec
        
        return (NoTangent(), Δcoord_params, NoTangent())
    end

    return y, SS_pullback
end

function ChainRulesCore.rrule(
    ::Type{Cartesian}, coord::StiefelScheifele, μ::Real; kwargs...
)
    # Forward evaluation
    y = Cartesian(coord, μ; kwargs...)
    
    function Cartesian_from_SS_pullback(Δy)
        # Extract parameters for the transformation
        coord_params = AstroCoords.params(coord)
        
        # Convert Tangent to vector if needed
        if Δy isa Tangent
            Δy_vec = collect(values(Δy))
        else
            Δy_vec = Δy
        end
        
        # Compute the Jacobian ∂(Cartesian)/∂(StiefelScheifele) analytically
        # Use Zygote to differentiate the underlying transformation function
        J = jacobian(x -> AstroCoords.StiefelScheifele2cart(x, μ; kwargs...), coord_params)
        
        # Apply chain rule: Δcoord = J^T * Δy
        Δcoord_params = J' * Δy_vec
        
        return (NoTangent(), Δcoord_params, NoTangent())
    end

    return y, Cartesian_from_SS_pullback
end

end
