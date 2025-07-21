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
        # We need the gradient w.r.t. each input argument
        grad_vals = collect(Δ)  # Convert to vector if needed
        
        # Each gradient component corresponds to one input argument
        return (NoTangent(), grad_vals...)
    end
    
    return result, astrocoord_pullback
end

end
