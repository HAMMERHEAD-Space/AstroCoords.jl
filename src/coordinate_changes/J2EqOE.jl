"""
    function koeM2IOE(u::AbstractVector{T}, μ::V) where {T<:Number, V<:Number}
        
Computes the Intermediate Orbit Elements from a Keplerian set.

# Arguments
-`u::AbstractVector{<:Number}`: The Keplerian state vector [a; e; i; Ω(RAAN); ω(AOP); M(Mean Anomaly)].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_IOE::SVector{6, <:Number}`: The Intermediate Orbit Element vector [I1; I2; I3; I4; I5; I6].
"""
function koe2IOE(u::AbstractVector{T}, μ::V) where {T<:Number, V<:Number}
    RT = promote_type(T, V)

    a, e, i, Ω, ω, f = u

    M = trueAnomaly2MeanAnomaly(f, e)

    I1 = a
    I2 = e * sin(M)
    I3 = e * cos(M)
    I4 = sin(i / 2.0) * sin(Ω)
    I5 = sin(i / 2.0) * cos(Ω)
    I6 = rem2pi(M + ω + Ω, RoundDown)

    return SVector{6,RT}(I1, I2, I3, I4, I5, I6)
end

"""
    function IOE2koeM(u::AbstractVector{T}, μ::V) where {T<:Number, V<:Number}

Computes the Keplerian orbital elements from a Intermediate Orbit Element set.

# Arguments
-`u::AbstractVector{<:Number}`: The Intermediate Orbit Element vector [I1; I2; I3; I4; I5; I6].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_koeM::SVector{6, <:Number}`: The Keplerian state vector [a; e; i; Ω(RAAN); ω(AOP); M(Mean Anomaly)].
"""
function IOE2koe(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    I1, I2, I3, I4, I5, I6 = u

    a = I1
    e = √(I2^2.0 + I3^2.0)
    i = 2.0 * asin(√(I4^2.0 + I5^2.0))
    Ω = atan(I4, I5)
    M = atan(I2, I3)
    ω = I6 - Ω - M

    f = meanAnomaly2TrueAnomaly(M, e)

    return SVector{6,RT}(a, e, i, Ω, ω, f)
end

"""
    function modEq2IOE(u::AbstractVector{T}, μ::V) where {T<:Number, V<:Number}

Computes the Intermediate Orbit Elements from a Modified Equinoctial set.

# Arguments
-`u::AbstractVector{<:Number}`: The Modified Equinoctial state vector [p; f; g; h; k; L].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_IOE::SVector{6, <:Number}`: The Intermediate Orbit Element vector [I1; I2; I3; I4; I5; I6].
"""
function modEqN2IOE(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    n, f, g, h, k, L = u

    I1 = ∛(μ / n^2.0)
    I2 = g * sin(L) - f * cos(L)
    I3 = f * sin(L) + g * cos(L)
    I4 = h / √(1.0 + h^2.0 + k^2.0)
    I5 = k / √(1.0 + h^2.0 + k^2.0)
    I6 = L

    return SVector{6,RT}(I1, I2, I3, I4, I5, I6)
end

"""
    function IOE2modEq(u::AbstractVector{T}, μ::V) where {T<:Number, V<:Number}

Computes the Modified Equinoctial orbital elements from a Intermediate Orbit Element set.

# Arguments
-`u::AbstractVector{<:Number}`: The Intermediate Orbit Element vector [I1; I2; I3; I4; I5; I6].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_modEq::SVector{6, <:Number}`: The Modified Equinoctial state vector [n; f; g; h; k; L].
"""
function IOE2modEqN(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    I1, I2, I3, I4, I5, I6 = u

    a = I1
    n = √(μ / (a^3.0))
    f = I3 * sin(I6) - I2 * cos(I6)
    g = I2 * sin(I6) + I3 * cos(I6)
    h = I4 / √(1.0 - I4^2.0 - I5^2.0)
    k = I5 / √(1.0 - I4^2.0 - I5^2.0)
    L = I6

    return SVector{6,RT}(n, f, g, h, k, L)
end

################################################################
# J2 Perturbed Equinoctial Orbital Elements
# Algorithm Steps
# These are broken out to test the individual steps
# See Appendix: Numerical Example
# Aristoff JM, Horwood JT, Alfriend KT (2021) On a set of J2 equinoctial orbital
# elements and their use for uncertainty propagation. Celestial Mechanics and
# Dynamical Astronomy 133(9):1–19
################################################################

@inline function _step1(u::AbstractVector{<:Number}, μ::Number)
    aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, fⱼ₂ = IOE2koe(u, μ)
    Lⱼ₂ = trueAnomaly2MeanAnomaly(fⱼ₂, eⱼ₂)

    return aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, fⱼ₂, Lⱼ₂
end

@inline function _step2(J2::Number, Req::Number, eⱼ₂::Number, aⱼ₂::Number, Iⱼ₂::Number)
    k = 0.5 * J2 * Req^2.0
    η = √(1.0 - eⱼ₂^2.0)
    γ = k / aⱼ₂^2.0
    γ′ = γ / η^4.0
    θ = cos(Iⱼ₂)

    return k, η, γ, γ′, θ
end

@inline function _step3(fⱼ₂::Number, eⱼ₂::Number)
    Eⱼ₂ = trueAnomaly2EccentricAnomaly(fⱼ₂, eⱼ₂)

    return Eⱼ₂
end

@inline function _step4(aⱼ₂::Number, eⱼ₂::Number, η::Number, Eⱼ₂::Number, Lⱼ₂::Number)
    rⱼ₂ = aⱼ₂ * (1.0 - eⱼ₂ * cos(Eⱼ₂))
    νⱼ₂ = atan(η * sin(Eⱼ₂), cos(Eⱼ₂) - eⱼ₂)
    #* Branch cut correction
    b_cut = νⱼ₂ - π
    Lⱼ₂ = Lⱼ₂ + ceil((b_cut - Lⱼ₂) / (2π))

    return rⱼ₂, νⱼ₂, Lⱼ₂
end

@inline function _step5(
    aⱼ₂::Number,
    γ::Number,
    γ′::Number,
    θ::Number,
    rⱼ₂::Number,
    η::Number,
    eⱼ₂::Number,
    gⱼ₂::Number,
    νⱼ₂::Number,
    Lⱼ₂::Number,
)
    a =
        aⱼ₂ * (
            1.0 +
            γ * (
                (-1.0 + 3.0 * θ^2.0) * ((aⱼ₂^3.0) / (rⱼ₂^3.0) - η^(-3.0)) +
                3.0 * (1.0 - θ^2.0) * (aⱼ₂^3.0) / (rⱼ₂^3.0) * cos(2.0 * gⱼ₂ + 2.0 * νⱼ₂)
            )
        )
    δh =
        -0.5 *
        γ′ *
        θ *
        (
            6.0 * (νⱼ₂ - Lⱼ₂ + eⱼ₂ * sin(νⱼ₂)) - 3.0 * sin(2.0 * gⱼ₂ + 2.0 * νⱼ₂) -
            3.0 * eⱼ₂ * sin(2.0 * gⱼ₂ + νⱼ₂) - eⱼ₂ * sin(2.0 * gⱼ₂ + 3.0 * νⱼ₂)
        )

    return a, δh
end

@inline function _step6(
    γ′::Number,
    θ::Number,
    νⱼ₂::Number,
    Lⱼ₂::Number,
    eⱼ₂::Number,
    gⱼ₂::Number,
    hⱼ₂::Number,
    δh::Number,
)
    Σlgh =
        0.25 *
        γ′ *
        (
            6.0 * (-1.0 + 5.0 * θ^2.0) * (νⱼ₂ - Lⱼ₂ + eⱼ₂ * sin(νⱼ₂)) +
            (3.0 - 5.0 * θ^2.0) * (
                3.0 * sin(2.0 * gⱼ₂ + 2.0 * νⱼ₂) +
                3.0 * eⱼ₂ * sin(2.0 * gⱼ₂ + νⱼ₂) +
                eⱼ₂ * sin(2.0 * gⱼ₂ + 3.0 * νⱼ₂)
            )
        ) +
        δh +
        Lⱼ₂ +
        gⱼ₂ +
        hⱼ₂

    return Σlgh
end

@inline function _step7(
    νⱼ₂::Number,
    eⱼ₂::Number,
    η::Number,
    aⱼ₂::Number,
    rⱼ₂::Number,
    γ::Number,
    γ′::Number,
    θ::Number,
    gⱼ₂::Number,
)
    v1 = 3.0 * cos(νⱼ₂) + 3.0 * eⱼ₂ * cos(νⱼ₂)^2.0 + eⱼ₂^2.0 * cos(νⱼ₂)^3.0
    v2 = η^(-6.0) * (eⱼ₂ * η + eⱼ₂ / (1.0 + η) + v1)
    v3 = η^(-6.0) * (eⱼ₂ + v1)
    v4 = ((aⱼ₂ * η) / rⱼ₂)^2.0 + aⱼ₂ / rⱼ₂

    δe =
        0.5 *
        η^2.0 *
        (
            γ * ((-1.0 + 3.0 * θ^2.0) * v2 + 3.0 * (1.0 - θ^2.0) * v3 * cos(2.0 * gⱼ₂ + 2.0 * νⱼ₂)) -
            γ′ * (1.0 - θ^2.0) * (3.0 * cos(2.0 * gⱼ₂ + νⱼ₂) + cos(2.0 * gⱼ₂ + 3.0 * νⱼ₂))
        )
    e″δL =
        -0.25 *
        η^3.0 *
        γ′ *
        (
            2.0 * (-1.0 + 3.0 * θ^2.0) * (v4 + 1.0) * sin(νⱼ₂) +
            3.0 *
            (1.0 - θ^2.0) *
            ((-v4 + 1.0) * sin(2.0 * gⱼ₂ + νⱼ₂) + (v4 + 1.0 / 3.0) * sin(2.0 * gⱼ₂ + 3.0 * νⱼ₂))
        )

    return v1, v2, v3, v4, δe, e″δL
end

@inline function _step8(
    γ′::Number, θ::Number, νⱼ₂::Number, eⱼ₂::Number, gⱼ₂::Number, δh::Number, Iⱼ₂::Number
)
    δI =
        0.5 *
        γ′ *
        θ *
        √(1.0 - θ^2.0) *
        (
            3.0 * cos(2.0 * gⱼ₂ + 2.0 * νⱼ₂) +
            3.0 * eⱼ₂ * cos(2.0 * gⱼ₂ + νⱼ₂) +
            eⱼ₂ * cos(2.0 * gⱼ₂ + 3.0 * νⱼ₂)
        )
    sin_half_I″_δh = (sin(Iⱼ₂) * δh) / (2.0 * cos(0.5 * Iⱼ₂))

    return δI, sin_half_I″_δh
end

@inline function _step9(
    a::Number,
    eⱼ₂::Number,
    δe::Number,
    e″δL::Number,
    Lⱼ₂::Number,
    Iⱼ₂::Number,
    δI::Number,
    sin_half_I″_δh::Number,
    hⱼ₂::Number,
    Σlgh::Number,
)
    I1 = a
    I2 = (eⱼ₂ + δe) * sin(Lⱼ₂) + e″δL * cos(Lⱼ₂)
    I3 = (eⱼ₂ + δe) * cos(Lⱼ₂) - e″δL * sin(Lⱼ₂)
    I4 = (sin(0.5 * Iⱼ₂) + cos(0.5 * Iⱼ₂) * 0.5 * δI) * sin(hⱼ₂) + sin_half_I″_δh * cos(hⱼ₂)
    I5 = (sin(0.5 * Iⱼ₂) + cos(0.5 * Iⱼ₂) * 0.5 * δI) * cos(hⱼ₂) - sin_half_I″_δh * sin(hⱼ₂)
    I6 = Σlgh

    return I1, I2, I3, I4, I5, I6
end

"""
    function J2IOE2IOE(u::AbstractVector{T}, μ::V) where {T<:Number, V<:Number}

Computes the Intermediate Orbit Elements from a J2 Perturbed Intermediate Orbit Element set.

# Arguments
-`u::AbstractVector{<:Number}`: The J2 Perturbed Intermediate Orbit Element vector [a; e; i; Ω(RAAN); ω(AOP); f(True Anomaly)].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_IOE::SVector{6, <:Number}`: The Intermediate Orbit Element vector [I1; I2; I3; I4; I5; I6].
"""
function J2IOE2IOE(
    u::AbstractVector{T}, μ::V; J2::W=1.0826261738522e-03, Req::Z=6.378137e+03
) where {T<:Number,V<:Number,W<:Number,Z<:Number}
    RT = promote_type(T, V, W, Z)

    aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, fⱼ₂, Lⱼ₂ = _step1(u, μ)
    _, η, γ, γ′, θ = _step2(J2, Req, eⱼ₂, aⱼ₂, Iⱼ₂)
    Eⱼ₂ = _step3(fⱼ₂, eⱼ₂)
    rⱼ₂, νⱼ₂, Lⱼ₂ = _step4(aⱼ₂, eⱼ₂, η, Eⱼ₂, Lⱼ₂)
    a, δh = _step5(aⱼ₂, γ, γ′, θ, rⱼ₂, η, eⱼ₂, gⱼ₂, νⱼ₂, Lⱼ₂)
    Σlgh = _step6(γ′, θ, νⱼ₂, Lⱼ₂, eⱼ₂, gⱼ₂, hⱼ₂, δh)
    _, _, _, _, δe, e″δL = _step7(νⱼ₂, eⱼ₂, η, aⱼ₂, rⱼ₂, γ, γ′, θ, gⱼ₂)
    δI, sin_half_I″_δh = _step8(γ′, θ, νⱼ₂, eⱼ₂, gⱼ₂, δh, Iⱼ₂)
    I1, I2, I3, I4, I5, I6 = _step9(
        a, eⱼ₂, δe, e″δL, Lⱼ₂, Iⱼ₂, δI, sin_half_I″_δh, hⱼ₂, Σlgh
    )

    return SVector{6,RT}(I1, I2, I3, I4, I5, I6)
end

"""
    function IOE2J2IOE(u::AbstractVector{T}, μ::V) where {T<:Number, V<:Number}

Computes the J2 Perturbed Intermediate Orbit Elements from a Intermediate Orbit Element set.

# Arguments
-`u::AbstractVector{<:Number}`: The Intermediate Orbit Element vector [I1; I2; I3; I4; I5; I6].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_J2IOR::SVector{6, <:Number}`: The J2 Perturbed Intermediate Orbit Element vector [a; e; i; Ω(RAAN); ω(AOP); f(True Anomaly)].
"""
function IOE2J2IOE(
    u::AbstractVector{T},
    μ::V;
    Req::W=6.378137e+03,
    tol::AbstractFloat=10.0 * eps(promote_type(T, V)),
    max_iter::Int=100,
) where {T<:Number,V<:Number,W<:Number}
    RT = promote_type(T, V, W)

    # Pre-allocate vectors for iteration
    scale = SVector{6,RT}(Req, 1.0, 1.0, 1.0, 1.0, 1.0)
    u_J2IOE_guess = SVector{6,RT}(u[1], u[2], u[3], u[4], u[5], u[6])
    u_IOE_guess = similar(u_J2IOE_guess)
    residual = similar(u_J2IOE_guess)

    # Initial guess setup
    u_J2IOE_guess = SVector{6,RT}(
        max(u_J2IOE_guess[1], 1e-2),
        clamp(u_J2IOE_guess[2], -1.0, 1.0),
        clamp(u_J2IOE_guess[3], -1.0, 1.0),
        clamp(u_J2IOE_guess[4], -1.0, 1.0),
        clamp(u_J2IOE_guess[5], -1.0, 1.0),
        u_J2IOE_guess[6],
    )

    # Iteration variables
    iter = 0
    error = typemax(RT)
    best_residual = error
    best_guess = u_J2IOE_guess

    while error > tol && iter < max_iter
        iter += 1

        # Compute IOE guess from J2IOE guess
        u_IOE_guess = J2IOE2IOE(u_J2IOE_guess, μ)
        u_IOE_guess = (
            u_IOE_guess[1] / scale[1],
            u_IOE_guess[2] / scale[2],
            u_IOE_guess[3] / scale[3],
            u_IOE_guess[4] / scale[4],
            u_IOE_guess[5] / scale[5],
            u_IOE_guess[6] / scale[6],
        )

        # Branch cut correction
        b_cut = u_J2IOE_guess[6] - π
        u_IOE_guess = SVector{6,RT}(
            u_IOE_guess[1],
            u_IOE_guess[2],
            u_IOE_guess[3],
            u_IOE_guess[4],
            u_IOE_guess[5],
            u_IOE_guess[6] + 2.0π * ceil((b_cut - u_IOE_guess[6]) / (2.0π)),
        )

        # Compute residual
        residual = (
            u[1] / scale[1] - u_IOE_guess[1],
            u[2] / scale[2] - u_IOE_guess[2],
            u[3] / scale[3] - u_IOE_guess[3],
            u[4] / scale[4] - u_IOE_guess[4],
            u[5] / scale[5] - u_IOE_guess[5],
            u[6] / scale[6] - u_IOE_guess[6],
        )

        # Update error and best guess if needed
        error = norm(residual)
        if error < best_residual
            best_residual = error
            best_guess = u_J2IOE_guess
        end

        # Update guess
        u_J2IOE_guess = SVector{6,RT}(
            u_J2IOE_guess[1] + residual[1] * Req,
            clamp(u_J2IOE_guess[2] + residual[2], -1.0, 1.0),
            clamp(u_J2IOE_guess[3] + residual[3], -1.0, 1.0),
            clamp(u_J2IOE_guess[4] + residual[4], -1.0, 1.0),
            clamp(u_J2IOE_guess[5] + residual[5], -1.0, 1.0),
            u_J2IOE_guess[6] + residual[6],
        )
    end

    return best_guess
end

"""
    function cart2J2EqOE(u::AbstractVector{<:Number}, μ::Number)

Computes the J2 Perturbed Equinoctial Orbit Elements from a Cartesian state vector.

# Arguments
-`u::AbstractVector{<:Number}`: The Cartesian state vector [x; y; z; ẋ; ẏ; ż].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_J2EqOE::SVector{6, <:Number}`: The J2 Perturbed Equinoctial Orbit Element vector [n; h; k; p; q; L].
"""
function cart2J2EqOE(u::AbstractVector{<:Number}, μ::Number)
    u_koe = cart2koe(u, μ)
    u_IOE = koe2IOE(u_koe, μ)
    u_J2IOE = IOE2J2IOE(u_IOE, μ)

    return IOE2modEqN(u_J2IOE, μ)
end

"""
    function J2EqOE2cart(u::AbstractVector{<:Number}, μ::Number)

Computes the Cartesian state vector from a J2 Perturbed Equinoctial Orbit Elements.

# Arguments
-`u::AbstractVector{<:Number}`: The J2 Perturbed Equinoctial Orbit Element vector [n; h; k; p; q; L].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_cart::SVector{6, <:Number}`: The Cartesian state vector [x; y; z; ẋ; ẏ; ż].
"""
function J2EqOE2cart(u::AbstractVector{<:Number}, μ::Number)
    u_J2IOE = modEqN2IOE(u, μ)
    u_IOE = J2IOE2IOE(u_J2IOE, μ)
    u_koe = IOE2koe(u_IOE, μ)

    return koe2cart(u_koe, μ)
end
