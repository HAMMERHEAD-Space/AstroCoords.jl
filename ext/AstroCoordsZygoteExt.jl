module AstroCoordsZygoteExt

using AstroCoords

using Zygote.ChainRulesCore: ChainRulesCore
import Zygote.ChainRulesCore: Tangent, NoTangent, ProjectTo

function ChainRulesCore.rrule(
    new_coord::Type{<:AstroCoords.AstroCoord}, coord::AbstractArray
)
    # Forward evaluation (Keplerian transformation)
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

end
