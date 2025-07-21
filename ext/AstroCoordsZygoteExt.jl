module AstroCoordsZygoteExt

using AstroCoords

using Zygote.ChainRulesCore: ChainRulesCore
import Zygote.ChainRulesCore: rrule, NoTangent

# Generic rule for any AstroCoord constructor
function ChainRulesCore.rrule(::Type{T}, args...) where {T<:AstroCoords.AstroCoord}
    # Forward evaluation
    result = T(args...)

    function astrocoord_pullback(Δ)
        # Convert gradient back to the parameter space
        # Δ should be compatible with params(result)
        grad_params = if Δ isa AbstractVector
            Δ
        else
            # Extract values if it's a structured gradient
            collect(Δ)
        end

        # Return NoTangent for the type and the gradient for the parameters
        return (NoTangent(), grad_params...)
    end

    return result, astrocoord_pullback
end

end
