"""
    function cart2koe(
        u::AbstractVector{T}, μ::V; equatorial_tol::Float64=1E-15, circular_tol::Float64=1E-15
    ) where {T<:Number,V<:Number}

Computes the Keplerian orbital elements from a cartesian set.

!!! note
    All angles are in radians.

# Arguments
-`u::AbstractVector{<:Number}`: The Cartesian state vector [x; y; z; ẋ; ẏ; ż].
-`μ::Number`: Standard graviational parameter of central body.

# Keyword Arguments
-`equatorial_tol::Float64`: The tolerance on what is considered an equatorial orbit (no inclination). [Default=1e-15]
-`circular_tol::Float64`: The tolerance on what is considered a circular orbit (no eccentricity). [Default=1e-15]

# Returns
-`u_koe::SVector{6, <:Number}``: Keplerian orbital element vector [a; e; i; Ω(RAAN); ω(AOP); f(True Anomaly)].
"""
function cart2koe(
    u::AbstractVector{T}, μ::V; equatorial_tol::Float64=1E-14, circular_tol::Float64=1E-14
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    x, y, z, ẋ, ẏ, ż = u

    r = SVector{3}(x, y, z)
    v = SVector{3}(ẋ, ẏ, ż)

    rmag = norm(r)
    vmag = norm(v)

    #* Angular Momentum
    h = cross(r, v)
    hmag = norm(h)

    #* Inclination
    k̂ = SVector{3,RT}(0.0, 0.0, 1.0)
    i = angle_between_vectors(h, k̂)

    #* Semi-Major Axis
    a = 1.0 / (2.0 / rmag - vmag^2 / μ)

    #* Eccentricity
    e = cross(v, h) / μ - r / rmag
    emag = norm(e)

    #* RAAN, AOP, True Anomaly
    if abs(i) < equatorial_tol
        if abs(emag) < circular_tol
            Ω = 0.0
            ω = 0.0
            f_raw = rem2pi(atan(y, x), RoundDown)
            # For circular equatorial, true longitude is well-defined but can wrap
            # If very close to 0 or 2π due to numerical errors, set to 0
            f = (abs(f_raw) < circular_tol || abs(f_raw - 2π) < circular_tol) ? 0.0 : f_raw
        else
            Ω = 0.0
            ω = rem2pi(atan(e[2], e[1]), RoundDown)
            f = rem2pi(atan(dot(h, cross(e, r)) / norm(h), dot(r, e)), RoundDown)
        end
    elseif abs(emag) < circular_tol
        n = cross(k̂, h)
        Ω = rem2pi(atan(n[2], n[1]), RoundDown)
        ω = 0.0
        f = rem2pi(atan(dot(r, cross(h, n)) / norm(h), dot(r, n)), RoundDown)
    else
        #* Eccentricity and True Anomaly
        if a > 0.0
            eSE = dot(r, v) / √(μ * a)
            eCE = (rmag * vmag^2) / μ - 1.0
            E = atan(eSE, eCE)
            f = rem2pi(eccentricAnomaly2TrueAnomaly(E, emag), RoundDown)
        else
            eSH = dot(r, v) / √(-μ * a)
            eCH = (rmag * vmag^2) / μ - 1.0
            F = log((eCH + eSH) / (eCH - eSH)) / 2.0
            f = rem2pi(
                2.0 * atan(√(1.0 + emag) * sinh(F / 2.0), √(emag - 1.0) * cosh(F / 2.0)),
                RoundDown,
            )
        end

        #* RAAN
        #* Line of Nodes
        n = cross(k̂, h)
        Ω = rem2pi(atan(n[2], n[1]), RoundDown)

        #* Argument of Perigee
        px = dot(r, n)
        py = dot(r, cross(h, n) / hmag)
        ω = rem2pi(atan(py, px) - f, RoundDown)
    end

    return SVector{6,RT}(a, emag, i, Ω, ω, f)
end

"""
    koe2cart(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Computes the Cartesian orbital elements from a Keplerian set.

!!! note
    All angles are in radians.

# Arguments
-`u::AbstractVector{<:Number}`: The Keplerian state vector [a; e; i; Ω(RAAN); ω(AOP); f(True Anomaly)].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_cart::Vector{6, <:Number}`: The Keplerian orbital element vector [x; y; z; ẋ; ẏ; ż].
"""
function koe2cart(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    a, e, i, Ω, ω, f = u

    rmag = (a * (1.0 - e^2) / (1.0 + e * cos(f)))

    θ = ω + f

    sΩ, cΩ = sincos(Ω)
    sθ, cθ = sincos(θ)
    sω, cω = sincos(ω)
    si, ci = sincos(i)

    x = rmag * (cΩ * cθ - sΩ * sθ * ci)
    y = rmag * (sΩ * cθ + cΩ * sθ * ci)
    z = rmag * (sθ * si)

    # Semi-latus rectum p = a(1-e²) is always positive for all orbit types
    # For elliptic: a > 0, e < 1, so a(1-e²) > 0
    # For hyperbolic: a < 0, e > 1, so a(1-e²) > 0 (both factors negative)
    # Angular momentum magnitude: h = √(μ*p)
    h = √(μ * a * (1.0 - e^2))

    ẋ = -μ / h * (cΩ * (sθ + e * sω) + sΩ * (cθ + e * cω) * ci)
    ẏ = -μ / h * (sΩ * (sθ + e * sω) - cΩ * (cθ + e * cω) * ci)
    ż = μ / h * (cθ + e * cω) * si

    return SVector{6,RT}(x, y, z, ẋ, ẏ, ż)
end

"""
    koe2USM7(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Converts Keplerian orbital elements into the Unified State Model set.
Van den Broeck, Michael. "An Approach to Generalizing Taylor Series Integration for Low-Thrust Trajectories." (2017).
https://repository.tudelft.nl/islandora/object/uuid%3A2567c152-ab56-4323-bcfa-b076343664f9

!!! note
    All angles are in radians.

# Arguments
-`u:AbstractVector{<:Number}`: The Keplerian state vector [a; e; i; Ω(RAAN); ω(AOP); f(True Anomaly)].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_USM::SVector{7, <:Number}`: The Unified State Model vector [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0].
"""
function koe2USM7(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    a, e, i, Ω, ω, f = u

    # Semi-latus rectum p = a(1-e²) is always positive for all orbit types
    # C = √(μ/p) where p is the semi-latus rectum
    C = √(μ / (a * (1.0 - e^2)))

    R = e * C
    Rf1 = -R * sin(Ω + ω)
    Rf2 = R * cos(Ω + ω)

    u = ω + f

    ϵO1 = sin(i / 2) * cos((Ω - u) / 2)
    ϵO2 = sin(i / 2) * sin((Ω - u) / 2)
    ϵO3 = cos(i / 2) * sin((Ω + u) / 2)
    η0 = cos(i / 2) * cos((Ω + u) / 2)

    return SVector{7,RT}(C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0)
end

"""
    USM72koe(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Converts Unified State Model elements into the Keplerian orbital element set.
Van den Broeck, Michael. "An Approach to Generalizing Taylor Series Integration for Low-Thrust Trajectories." (2017).
https://repository.tudelft.nl/islandora/object/uuid%3A2567c152-ab56-4323-bcfa-b076343664f9

!!! note
    All angles are in radians.

# Arguments
-`u::AbstractVector{<:Number}`: The Unified State Model vector [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_koe:SVector{6, <:Number}`: The Keplerian State vector [a; e; i; Ω(RAAN); ω(AOP); f(True Anomaly)].
"""
function USM72koe(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u

    sinλ = (2.0 * ϵO3 * η0) / (ϵO3^2 + η0^2)
    cosλ = (η0^2 - ϵO3^2) / (ϵO3^2 + η0^2)
    λ = atan(sinλ, cosλ)

    ve1 = Rf1 * cosλ + Rf2 * sinλ
    ve2 = C - Rf1 * sinλ + Rf2 * cosλ

    R = √(Rf1^2 + Rf2^2)
    e = R / C

    a = μ / (2 * C * ve2 - (ve1^2 + ve2^2))
    i = acos(1.0 - 2.0 * (ϵO1^2 + ϵO2^2))

    # Check for equatorial orbit (i ≈ 0, so ϵO1 ≈ 0 and ϵO2 ≈ 0)
    if (abs(ϵO1) < 1e-10 && abs(ϵO2) < 1e-10) || (abs(ϵO3) < 1e-10 && abs(η0) < 1e-10)
        Ω = 0.0
    else
        Ω = atan(
            (ϵO1 * ϵO3 + ϵO2 * η0) / √((ϵO1^2 + ϵO2^2) * (η0^2 + ϵO3^2)),
            (ϵO1 * η0 - ϵO2 * ϵO3) / √((ϵO1^2 + ϵO2^2) * (η0^2 + ϵO3^2)),
        )
        while Ω < 0.0
            Ω = Ω + 2 * π
        end
    end

    # Check for circular orbit (e ≈ 0, so R ≈ 0)
    if abs(R) < 1e-10
        ω = 0.0
        f = λ - Ω
        while f < 0.0
            f = f + 2 * π
        end
    else
        f = atan(ve1 / R, (ve2 - C) / R)
        while f < 0.0
            f = f + 2 * π
        end
        ω = λ - Ω - f
        while ω < 0.0
            ω = ω + 2 * π
        end
    end

    return SVector{6,RT}(a, e, i, Ω, ω, f)
end

"""
    USM72USM6(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Converts USM with quaternions to USM with Modified Rodrigue Parameters.

# Arguments
-`u::AbstractVector{<:Number}`: The Unified State Model vector [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_USM6::SVector{6, <:Number}`: The USM6 State vector [C; Rf1; Rf2; σ1; σ2; σ3].
"""
function USM72USM6(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u

    EPs = SVector{4}(η0, ϵO1, ϵO2, ϵO3)

    σ = EP2MRP(EPs)

    return SVector{6,RT}(C, Rf1, Rf2, σ[1], σ[2], σ[3])
end

"""
    USM62USM7(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Converts USM with Modified Rodrigue Parameters to USM with quaternions.

# Arguments
-`u::AbstractVector{<:Number}`: The USM6 vector [C; Rf1; Rf2; σ1; σ2; σ3].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_USM::SVector{6, <:Number}`: The Unified State Model vector [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0].
"""
function USM62USM7(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    C, Rf1, Rf2, σ1, σ2, σ3 = u

    MRPs = SVector{3}(σ1, σ2, σ3)

    η0, ϵO1, ϵO2, ϵO3 = MRP2EP(MRPs)

    return SVector{7,RT}(C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0)
end

"""
    USM72USMEM(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Converts USM with quaternions to USM with exponential mapping.

# Arguments
-`u::AbstractVector{<:Number}`: The Unified State Model vector [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_USMEM::SVector{6, <:Number}`: The USMEM State Vector [C; Rf1; Rf2; a1; a2; a3].
"""
function USM72USMEM(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u

    Φ = 2.0 * acos(clamp(η0, -1.0, 1.0))

    # Handle singularity for zero rotation (equatorial orbit)
    if abs(Φ) < 1e-10
        # For small angles, ϵ ≈ a/2, so a ≈ 2ϵ and em = Φ*a ≈ 0
        return SVector{6,RT}(C, Rf1, Rf2, zero(RT), zero(RT), zero(RT))
    end

    denom = sin(Φ / 2.0)
    a = SVector{3}(ϵO1 / denom, ϵO2 / denom, ϵO3 / denom)

    em = Φ * a

    return SVector{6,RT}(C, Rf1, Rf2, em[1], em[2], em[3])
end

"""
    USMEM2USM7(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Converts USM with exponential mapping to USM with quaternions.

# Arguments
-`u::AbstractVector{<:Number}`: The USMEM vector [C; Rf1; Rf2; a1; a2; a3].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_USM::SVector{7, <:Number}`: The Unified State Model vector [C; Rf1; Rf2; ϵO1; ϵO2; ϵO3; η0].
"""
function USMEM2USM7(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    C, Rf1, Rf2, a1, a2, a3 = u

    a = SVector{3}(a1, a2, a3)
    Φ = norm(a)

    # Handle singularity for zero rotation (equatorial orbit)
    if abs(Φ) < 1e-10
        # For small Φ, sin(Φ/2)/Φ → 1/2, so ϵ → a/2
        ϵ = 0.5 * a
        η0 = one(RT)  # cos(0) = 1
    else
        ϵ = sin(Φ / 2.0) / Φ * a
        η0 = cos(Φ / 2.0)
    end

    return SVector{7,RT}(C, Rf1, Rf2, ϵ[1], ϵ[2], ϵ[3], η0)
end

"""
    koe2ModEq(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Converts Keplerian elements into the Modified Equinoctial elements.

!!! note
    All angles are in radians.

# Arguments
-`u:AbstractVector{<:Number}`: The Keplerian State vector [a; e; i; Ω(RAAN); ω(AOP); ν(True Anomaly)].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_ModEq::SVector{6, <:Number}`: The Modified Equinoctial state vector [p; f; g; h; k; l].
"""
function koe2ModEq(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    a, e, i, Ω, ω, ν = u

    p = a * (1 - e^2)
    f = e * cos(ω + Ω)
    g = e * sin(ω + Ω)
    h = tan(i / 2) * cos(Ω)
    k = tan(i / 2) * sin(Ω)
    L = Ω + ω + ν

    return SVector{6,RT}(p, f, g, h, k, L)
end

"""
    ModEq2koe(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Converts Modified Equinoctial elements into the Keplerian elements.

!!! note
    All angles are in radians.

# Arguments
-`u:AbstractVector{<:Number}`: The Modified Equinoctial state vector [p; f; g; h; k; l].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_koe::SVector{6, <:Number}`: The Keplerian state vector [a; e; i; Ω(RAAN); ω(AOP); ν(True Anomaly)].
"""
function ModEq2koe(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    p, f, g, h, k, L = u

    a = p / (1 - f^2 - g^2)
    e = √(f^2 + g^2)
    i = atan(2 * √(h^2 + k^2), 1 - h^2 - k^2)

    # Handle equatorial orbit singularity (i ≈ 0, h ≈ 0, k ≈ 0)
    # For equatorial orbits, Ω and ω are not individually defined
    # but their sum (longitude of periapsis) is well-defined
    equatorial = abs(h) < 1e-10 && abs(k) < 1e-10
    circular = abs(e) < 1e-10

    if equatorial && circular
        # Circular equatorial: only true longitude is defined
        Ω = zero(RT)
        ω = zero(RT)
        ν = L
    elseif equatorial
        # Eccentric equatorial: Ω undefined, use longitude of periapsis
        Ω = zero(RT)
        ω = atan(g, f)
        ν = L - ω
    elseif circular
        # Circular inclined: ω undefined, use argument of latitude
        Ω = atan(k, h)
        ω = zero(RT)
        ν = L - Ω
    else
        # General case: all elements well-defined
        Ω = atan(k, h)
        ω = atan(g * h - f * k, f * h + g * k)
        ν = L - Ω - ω
    end

    return SVector{6,RT}(a, e, i, Ω, ω, ν)
end

"""
    ModEq2ModEqN(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Converts Modified Equinoctial elements to Modified Equinoctial with mean motion.

!!! note
    All angles are in radians.

# Arguments
-`u:AbstractVector{<:Number}`: The Modified Equinoctial state vector [p; f; g; h; k; L].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_ModEqN::SVector{6, <:Number}`: The Modified Equinoctial (mean motion) state vector [η; f; g; h; k; L].
"""
function ModEq2ModEqN(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    p, f, g, h, k, L = u

    # Compute eccentricity squared
    e_sq = f^2 + g^2
    
    # Compute semi-major axis: a = p / (1 - e²)
    a = p / (1 - e_sq)
    
    # Compute mean motion: 
    # For elliptic orbits (a > 0, e < 1): n = √(μ/a³)
    # For hyperbolic orbits (a < 0, e > 1): n = √(μ/|a|³) = √(-μ/a³)
    if a > 0
        η = √(μ / a^3)
    else
        η = √(-μ / a^3)
    end

    return SVector{6,RT}(η, f, g, h, k, L)
end

"""
    ModEqN2ModEq(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Converts Modified Equinoctial with mean motion to Modified Equinoctial elements.

!!! note
    All angles are in radians.

# Arguments
-`u:AbstractVector{<:Number}`: The Modified Equinoctial (mean motion) state vector [η; f; g; h; k; L].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_ModEq::SVector{6, <:Number}`: The Modified Equinoctial state vector [p; f; g; h; k; L].
"""
function ModEqN2ModEq(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    η, f, g, h, k, L = u

    # Compute eccentricity squared and eccentricity
    e_sq = f^2 + g^2
    e = √(e_sq)
    
    # Compute semi-major axis from mean motion:
    # For elliptic orbits (e < 1): a = (μ/n²)^(1/3) > 0
    # For hyperbolic orbits (e > 1): a = -(μ/n²)^(1/3) < 0
    a_mag = ∛(μ / η^2)
    if e < 1.0
        # Elliptic orbit
        a = a_mag
    else
        # Hyperbolic orbit
        a = -a_mag
    end
    
    # Compute semi-parameter: p = a(1 - e²)
    # For hyperbolic: a < 0, e² > 1, so (1 - e²) < 0, making p > 0
    p = a * (1 - e_sq)

    return SVector{6,RT}(p, f, g, h, k, L)
end

"""
    cart2Mil(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Converts Cartesian state vector into the Milankovich state vector.

# Arguments
-`u::AbstractVector{<:Number}`: The Cartesian state vector [x; y; z; ẋ; ẏ; ż].
-`μ::Number`: Standard graviational parameter of central body.

# Keyword Arguments
-`equatorial_tol::Float64`: The tolerance on what is considered an equatorial orbit (no inclination). [Default=1e-15]
-`circular_tol::Float64`: The tolerance on what is considered a circular orbit (no eccentricity). [Default=1e-15]


# Returns
-`u_Mil::SVector{7, <:Number}`: The Milankovich state vector [H; e; L]. 
"""
function cart2Mil(
    u::AbstractVector{T}, μ::V; equatorial_tol::Float64=1E-15, circular_tol::Float64=1E-15
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    r = SVector{3}(u[1], u[2], u[3])
    v = SVector{3}(u[4], u[5], u[6])

    H = cross(r, v)
    e = cross(v / μ, H) - r / norm(r)

    _, _, _, Ω, ω, f = cart2koe(
        u, μ; equatorial_tol=equatorial_tol, circular_tol=circular_tol
    )
    L = Ω + ω + f

    return SVector{7,RT}(H[1], H[2], H[3], e[1], e[2], e[3], L)
end

"""
   Mil2cart(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Converts Milankovich state vector into the Cartesian state vector.

# Arguments
-`u::AbstractVector{<:Number}`: The Milankovich state vector [H; e; L].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_cart::SVector{6, <:Number}`: The Cartesian state vector [x; y; z; ẋ; ẏ; ż].
"""
function Mil2cart(
    u::AbstractVector{T}, μ::V; equatorial_tol::Float64=1E-15, circular_tol::Float64=1E-15
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    H = SVector{3}(u[1], u[2], u[3])
    e = SVector{3}(u[4], u[5], u[6])
    L = u[7]

    ẑ = SVector{3}(0.0, 0.0, 1.0)

    i = angle_between_vectors(H, ẑ)

    hmag = √(sum(abs2.(H)))
    emag = √(sum(abs2.(e)))

    if abs(i) < equatorial_tol
        if abs(emag) < circular_tol
            Ω = 0.0
            ω = 0.0
        else
            Ω = 0
            ω = rem2pi(atan(e[2], e[1]), RoundDown)
        end
    elseif abs(emag) < circular_tol
        n = cross(ẑ, H ./ hmag)
        Ω = rem2pi(atan(n[2], n[1]), RoundDown)
        ω = 0.0
    else
        n = cross(ẑ, H ./ hmag)
        Ω = rem2pi(atan(n[2], n[1]), RoundDown)
        ω = angle_between_vectors(n, e)
    end

    f = L - Ω - ω

    a = (hmag^2) / (μ * (1 - emag^2))
    rmag = (a * (1.0 - emag^2) / (1.0 + emag * cos(f)))

    θ = ω + f

    sΩ, cΩ = sincos(Ω)
    sθ, cθ = sincos(θ)
    sω, cω = sincos(ω)
    si, ci = sincos(i)

    x = rmag * (cΩ * cθ - sΩ * sθ * ci)
    y = rmag * (sΩ * cθ + cΩ * sθ * ci)
    z = rmag * (sθ * si)

    ẋ = -μ / hmag * (cΩ * (sθ + emag * sω) + sΩ * (cθ + emag * cω) * ci)
    ẏ = -μ / hmag * (sΩ * (sθ + emag * sω) - cΩ * (cθ + emag * cω) * ci)
    ż = μ / hmag * (cθ + emag * cω) * si

    return SVector{6,RT}(x, y, z, ẋ, ẏ, ż)
end

"""
   cart2cylind(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Computes the cylindrical orbital elements from a Cartesian set.

!!! note
    All angles are in radians.

# Arguments
-`u::AbstractVector{<:Number}`: Cartesian state vector [x; y; z; ẋ; ẏ; ż].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_cylind::SVector{6, <:Number}``: The cylindrical orbital element vector [r; θ; z; ṙ; θdot; ż].
"""
function cart2cylind(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    x, y, z, ẋ, ẏ, ż = u

    ρ_vec = SVector{2}(x, y)

    ρ = norm(ρ_vec)
    θ = atan(y, x)

    ρdot = (x * ẋ + y * ẏ) / ρ
    θdot = (x * ẏ - y * ẋ) / ρ

    return SVector{6,RT}(ρ, θ, z, ρdot, θdot, ż)
end

"""
    cylind2cart(u::AbstractVector{T}, μ::Number) where {T<:Number}

Computes the Cartesian orbital elements from a cylindrical set.

!!! note
    All angles are in radians.

# Arguments
-`u::AbstractVector{<:Number}`: The cylindrical state vector [r; θ; z; ṙ; θdot; ż].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_cart::SVector{6, <:Number}`: The Cartesian orbital element vector [x; y; z; ẋ; ẏ; ż].
"""
function cylind2cart(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    ρ, θ, z, ρdot, θdot, ż = u

    x = ρ * cos(θ)
    y = ρ * sin(θ)

    ẋ = ρdot * cos(θ) - θdot * sin(θ)
    ẏ = ρdot * sin(θ) + θdot * cos(θ)

    return SVector{6,RT}(x, y, z, ẋ, ẏ, ż)
end

"""
    cart2sphere(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Computes the spherical orbital elements from a spherical set.

!!! note
    All angles are in radians.

# Arguments
-`u::AbstractVector{<:Number}`: The Cartesian state vector [x; y; z; ẋ; ẏ; ż].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-'u_sphere::SVector{6, <:Number}': Spherical Orbital Element Vector [r; θ; ϕ; ṙ; θdot; ϕdot]
"""
function cart2sphere(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    x, y, z, ẋ, ẏ, ż = u

    r_vec = SVector{3}(x, y, z)

    r = norm(r_vec)
    θ = atan(y, x)
    ϕ = acos(z / r)

    ṙ = (x * ẋ + y * ẏ + z * ż) / r
    θdot = (x * ẏ - ẋ * y) / (x^2 + y^2)
    ϕdot = (z * (x * ẋ + y * ẏ) - ż * (x^2 + y^2)) / (√(x^2.0 + y^2.0) * r^2)

    return SVector{6,RT}(r, θ, ϕ, ṙ, θdot, ϕdot)
end

"""
    sphere2cart(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Computes the Cartesian orbital elements from a spherical set.

!!! note
    All angles are in radians.

# Arguments
-`u::AbstractVector{<:Number}`: The spherical orbital element vector [r; θ; ϕ; ṙ; θdot; ϕdot].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-'u_cart::SVector{6, <:Number}': The Cartesian orbital element vector [x; y; z; ẋ; ẏ; ż].
"""
function sphere2cart(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    r, θ, ϕ, ṙ, θdot, ϕdot = u

    x = r * cos(θ) * sin(ϕ)
    y = r * sin(θ) * sin(ϕ)
    z = r * cos(ϕ)

    ẋ = ṙ * cos(θ) * sin(ϕ) - r * θdot * sin(θ) * sin(ϕ) + r * ϕdot * cos(θ) * cos(ϕ)
    ẏ = ṙ * sin(θ) * sin(ϕ) + r * θdot * cos(θ) * sin(ϕ) + r * ϕdot * sin(θ) * cos(ϕ)
    ż = ṙ * cos(ϕ) - r * ϕdot * sin(ϕ)

    return SVector{6,RT}(x, y, z, ẋ, ẏ, ż)
end

"""
    cart2delaunay(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}

Computes the Delaunay orbital elements from a Cartesian set.
Laskar, Jacques. "Andoyer construction for Hill and Delaunay variables." Celestial Mechanics and Dynamical Astronomy 128.4 (2017): 475-482.

!!! note
    All angles are in radians.

# Arguments
-`u::AbstractVector{<:Number}`: The Cartesian orbital element vector [x; y; z; ẋ; ẏ; ż].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-'u_cart::SVector{6, <:Number}': Delaunay Orbital Element Vector [L; G; H; M; ω; Ω]
"""
function cart2delaunay(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    RT = promote_type(T, V)
    a, e, _, Ω, ω, f = cart2koe(u, μ)
    M = trueAnomaly2MeanAnomaly(f, e)

    r = SVector{3}(u[1], u[2], u[3])
    v = SVector{3}(u[4], u[5], u[6])

    h = cross(r, v)

    # For elliptic orbits (a > 0, e < 1): L = √(μ*a), e = √(1 - (G/L)²)
    # For hyperbolic orbits (a < 0, e > 1): L = √(-μ*a), e = √(1 + (G/L)²)
    if a > 0
        L = √(μ * a)
    else
        L = √(-μ * a)
    end
    G = norm(h)
    H = h[3]

    return SVector{6,RT}(L, G, H, M, ω, Ω)
end

"""
    delaunay2cart(u::AbstractVector{T}, μ::V) where {T<:Number,V<:Number}
    
Computes the Cartesian orbital elements from a Delaunay set.
Laskar, Jacques. "Andoyer construction for Hill and Delaunay variables." Celestial Mechanics and Dynamical Astronomy 128.4 (2017): 475-482.

!!! note
    All angles are in radians.


# Argument
-`u::AbstractVector{<:Number}`: The Delaunay orbital element vector [L; G; H; M; ω; Ω].
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`u_cart::SVector{6, <:Number}``: The cartesian orbital element vector [x; y; z; ẋ; ẏ; ż].
"""
function delaunay2cart(
    u::AbstractVector{T}, μ::V; circular_tol::Float64=1e-14
) where {T<:Number,V<:Number}
    RT = promote_type(T, V)

    L, G, H, M, ω, Ω = u

    # Check if orbit is hyperbolic: for hyperbolic orbits, G/L > 1
    # Elliptic/Circular: L = √(μ*a), e = √(1 - (G/L)²), a > 0
    # Hyperbolic: L = √(-μ*a), e = √(1 + (G/L)²), a < 0
    # Note: For circular orbits, G/L = 1 exactly
    G_over_L_sq = (G / L)^2

    # Check if nearly circular (to handle floating-point errors)
    if abs(√(G_over_L_sq) - 1.0) < circular_tol
        # Circular orbit: treat as exactly circular
        a = L^2 / μ
        e = 0.0
    elseif √(G_over_L_sq) < 1.0
        # Elliptic orbit
        a = L^2 / μ
        e = √(1.0 - G_over_L_sq)
    else
        # Hyperbolic orbit
        a = -L^2 / μ
        e = √(1.0 + G_over_L_sq)
    end

    i = acos(H / G)

    f = meanAnomaly2TrueAnomaly(M, e)

    u_koe = SVector{6,RT}(a, e, i, Ω, ω, f)

    return koe2cart(u_koe, μ)
end
