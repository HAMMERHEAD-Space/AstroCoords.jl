module AstroCoordsZygoteExt

using AstroCoords

using Zygote.ChainRulesCore: ChainRulesCore
import Zygote.ChainRulesCore: rrule, NoTangent

# Generic rule for any AstroCoord constructor
function ChainRulesCore.rrule(::Type{T}, args...) where {T<:AstroCoords.AstroCoord}
    # Forward evaluation
    result = T(args...)

    function astrocoord_pullback(Δ)
        # The gradient Δ is w.r.t. the output AstroCoord object
        # We need to extract the gradient for each input parameter
        if Δ isa AbstractVector
            # Direct vector gradient
            return (NoTangent(), Δ...)
        else
            # Try to extract from structured gradient (Tangent, etc.)
            # Convert to params vector first to get proper structure
            param_grad = params(Δ)
            return (NoTangent(), param_grad...)
        end
    end

    return result, astrocoord_pullback
end

end
