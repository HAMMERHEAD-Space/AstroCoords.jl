@doc raw"""
    ModEq{T} <: CoordinateSet{T}

Cartesian coordinates.
"""
struct ModEq{T} <: CoordinateSet{T}
    g::SVector{6,T}
end

function (::Type{J})(g::StaticArray{Tuple{N},T,1}) where {N,T,J<:ModEq}
    return J(SVector{6,T}(g))
end

ModEq(X::AbstractVector{T}) where {T} = ModEq(SVector{6,T}(X))

params(c::ModEq) = c.g
Base.eltype(::Type{ModEq{T}}) where {T} = T

function Base.show(io::IO, c::ModEq{T}) where {T}
    return print(io, "ModEq{$T}($(c.g))")
end


@doc raw"""
    J2EqOE{T} <: CoordinateSet{T}

J2-equinoctial orbital elements.
"""
struct J2EqOE{T} <: CoordinateSet{T}
    g::SVector{6,T}
end

function (::Type{J})(g::StaticArray{Tuple{N},T,1}) where {N,T,J<:J2EqOE}
    return J(SVector{6,T}(g))
end

J2EqOE(X::AbstractVector{T}) where {T} = J2EqOE(SVector{6,T}(X))

params(c::J2EqOE) = c.g
Base.eltype(::Type{J2EqOE{T}}) where {T} = T

function Base.show(io::IO, c::J2EqOE{T}) where {T}
    return print(io, "J2EqOE{$T}($(c.g))")
end


function Cartesian(c::ModEq, μ)
    n, h, k, p, q, L = c.g
    return modEqN2Cartesian(n, h, k, p, q, L, μ)
end

function ModEq(c::Cartesian, μ)
    r = c.g[SVector(1, 2, 3)]
    v = c.g[SVector(4, 5, 6)]
    return cartesian2ModEqN(r, v, μ)
end

function Cartesian(c::J2EqOE, μ, R=6378.137, J2=0.0010826261738522)
    # Convert to ModEqN coordinates
    I = modEqIOE2N(c.g, μ, R, J2)

    # Now use the ModEqN to Cartesian conversion
    n, h, k, p, q, L = I
    return modEqN2Cartesian(n, h, k, p, q, L, μ)
end

function J2EqOE(c::Cartesian, μ, R=6378.137, J2=0.0010826261738522)
    # Convert Cartesian to ModEqN coordinates
    r = c.g[SVector(1, 2, 3)]
    v = c.g[SVector(4, 5, 6)]
    modeqN = cartesian2ModEqN(r, v, μ)
    I = params(modeqN)

    # Convert from ModEqN to J2EqOE
    O = modEqN2IOE(I, μ, R, J2)
    return J2EqOE(O)
end

function Keplerian(c::ModEq, μ)
    n, h, k, p, q, L = c.g
    return modEqN2Keplerian(n, h, k, p, q, L, μ)
end

function ModEq(c::Keplerian, μ)
    a, e, i, Ω, ω, ν = c.g
    return keplerian2ModEqN(a, e, i, Ω, ω, ν)
end

function Keplerian(c::J2EqOE, μ, R=6378.137, J2=0.0010826261738522)
    # Convert J2EqOE to ModEqN
    I = modEqIOE2N(c.g, μ, R, J2)

    # Convert ModEqN to Keplerian
    n, h, k, p, q, L = I
    return modEqN2Keplerian(n, h, k, p, q, L, μ)
end

function J2EqOE(c::Keplerian, μ, R=6378.137, J2=0.0010826261738522)
    # Convert Keplerian to ModEqN
    a, e, i, Ω, ω, ν = c.g
    modeqN = keplerian2ModEqN(a, e, i, Ω, ω, ν)
    I = params(modeqN)

    # Convert ModEqN to J2EqOE
    O = modEqN2IOE(I, μ, R, J2)
    return J2EqOE(O)
end

function Delaunay(c::ModEq, μ)
    n, h, k, p, q, L = c.g
    return modEqN2Delaunay(n, h, k, p, q, L, μ)
end

function ModEq(c::Delaunay, μ)
    L, G, H, l, g, h = c.g
    return delaunay2ModEqN(L, G, H, l, g, h)
end

function Delaunay(c::J2EqOE, μ, R=6378.137, J2=0.0010826261738522)
    # Convert J2EqOE to ModEqN
    I = modEqIOE2N(c.g, μ, R, J2)

    # Convert ModEqN to Delaunay
    n, h, k, p, q, L = I
    return modEqN2Delaunay(n, h, k, p, q, L, μ)
end

function J2EqOE(c::Delaunay, μ, R=6378.137, J2=0.0010826261738522)
    # Convert Delaunay to ModEqN
    L, G, H, l, g, h = c.g
    modeqN = delaunay2ModEqN(L, G, H, l, g, h)
    I = params(modeqN)

    # Convert ModEqN to J2EqOE
    O = modEqN2IOE(I, μ, R, J2)
    return J2EqOE(O)
end

# Helper conversion functions
function modEqN2IOE(u, μ, R=6378.137, J2=0.0010826261738522)
    n, h, k, p, q, L = u

    # Compute the osculating Keplerian elements from the ModEqN coordinates
    kep = modEqN2Keplerian(n, h, k, p, q, L, μ)
    a, e, i, Ω, ω, ν = kep.g

    # Compute the J2-perturbed elements through the inverse transformation
    aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, _, Lⱼ₂ = _step1([a, h, k, p, q, L], μ)
    k, η, γ, γ′, θ = _step2(J2, R, eⱼ₂, aⱼ₂, Iⱼ₂)

    f = ν + ω
    fⱼ₂ = f - θ * sin(f + gⱼ₂)
    Eⱼ₂ = _step3(fⱼ₂, eⱼ₂)
    rⱼ₂, νⱼ₂, _ = _step4(aⱼ₂, eⱼ₂, η, Eⱼ₂, Lⱼ₂)
    aᵢₒₑ, δh = _step5(aⱼ₂, γ, γ′, θ, rⱼ₂, η, eⱼ₂, gⱼ₂, νⱼ₂, Lⱼ₂)
    Σlgh = _step6(γ′, θ, νⱼ₂, Lⱼ₂, eⱼ₂, gⱼ₂, hⱼ₂, δh)
    _, _, _, _, δe, e″δL = _step7(νⱼ₂, eⱼ₂, η, aⱼ₂, rⱼ₂, γ, γ′, θ, gⱼ₂)
    δI, sin_half_I″_δh = _step8(γ′, θ, νⱼ₂, eⱼ₂, gⱼ₂, δh, Iⱼ₂)
    I1, I2, I3, I4, I5, I6 = _step9(
        aᵢₒₑ, eⱼ₂, δe, e″δL, Lⱼ₂, Iⱼ₂, δI, sin_half_I″_δh, hⱼ₂, Σlgh
    )

    return SVector{6}(I1, I2, I3, I4, I5, I6)
end

function modEqIOE2N(u, μ, R=6378.137, J2=0.0010826261738522)
    a, h, k, p, q, L = u

    # Step 1: Convert IOE to J2-perturbed elements
    aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, _, Lⱼ₂ = _step1(u, μ)

    # Step 2: Compute auxiliary quantities
    kval, η, γ, γ′, θ = _step2(J2, R, eⱼ₂, aⱼ₂, Iⱼ₂)

    # Step 3-4: Compute true anomaly and radius
    fⱼ₂ = mod2pi(Lⱼ₂ - gⱼ₂ - hⱼ₂)
    Eⱼ₂ = _step3(fⱼ₂, eⱼ₂)
    rⱼ₂, νⱼ₂, uⱼ₂ = _step4(aⱼ₂, eⱼ₂, η, Eⱼ₂, Lⱼ₂)

    # Step 10: Compute the osculating elements (inverse of steps 5-9)
    n_osc, h_osc, k_osc, p_osc, q_osc, L_osc = _step10(
        a, h, k, p, q, L, aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, fⱼ₂, Lⱼ₂, γ, γ′, θ, rⱼ₂, η, νⱼ₂
    )

    return SVector{6}(n_osc, h_osc, k_osc, p_osc, q_osc, L_osc)
end

function _step1(u, μ)
    a, h, k, p, q, L = u

    n = sqrt(μ / a^3)

    # Compute eccentricity and inclination
    e = sqrt(h^2 + k^2)
    I = 2 * atan(sqrt(p^2 + q^2))

    # Compute auxiliary angles
    g = atan(k, h)
    h_angle = atan(q, p)
    f = mod2pi(L - g - h_angle)

    return a, e, I, h_angle, g, f, L
end

function _step2(J2, R, e, a, I)
    k = 3 * J2 * R^2 / 4
    η = sqrt(1 - e^2)
    γ = k / (a^2 * η^4)
    γ′ = k / (a^2 * η^3)
    θ = (1 - 3 * cos(I)^2 / 2) * γ

    return k, η, γ, γ′, θ
end

function _step3(f, e)
    sinf = sin(f)
    cosf = cos(f)

    E = atan(sqrt(1 - e^2) * sinf, e + cosf)

    return E
end

function _step4(a, e, η, E, L)
    r = a * (1 - e * cos(E))

    # Compute true anomaly from eccentric anomaly
    sinE = sin(E)
    cosE = cos(E)
    ν = atan(sqrt(1 - e^2) * sinE, cosE - e)

    u = mod2pi(ν + (L - ν))

    return r, ν, u
end

function _step5(a, γ, γ′, θ, r, η, e, g, ν, L)
    # Compute the semi-major axis correction
    δa = -2 * γ * r * (1 + (1 + η) * sin(ν + g)^2)
    a_osc = a + δa

    # Compute the ascending node correction
    δh = γ′ * sin(ν + g) * cos(ν + g) / tan(acos(sqrt(1 - (sin(L - ν - g))^2)))

    return a_osc, δh
end

function _step6(γ′, θ, ν, L, e, g, h, δh)
    Σlgh = L + θ * (1 + 1 / (1 - e^2)) * sin(ν + g)^2 - γ′ * sin(ν + g) * cos(ν + g) / tan(acos(sqrt(1 - (sin(L - ν - g))^2))) + δh

    return Σlgh
end

function _step7(ν, e, η, a, r, γ, γ′, θ, g)
    v1 = cos(ν + g)
    v2 = cos(2 * ν + 2 * g)
    v3 = cos(ν + g)^2
    v4 = 2 * ν + 2 * g

    δe = (2 * γ / (a * η)) * (1 + η) * sin(ν + g) * v1
    e″δL = -θ * ((2 + e^2) / e) * v2 - γ′ * v1 / (e * tan(acos(sqrt(1 - (sin(ν + g))^2))))

    return v1, v2, v3, v4, δe, e″δL
end

function _step8(γ′, θ, ν, e, g, δh, I)
    δI = -2 * γ′ * cos(ν + g) * sin(ν + g) / sin(I)

    sin_half_I″_δh = γ′ * cos(ν + g) * sin(ν + g) / tan(I / 2)

    return δI, sin_half_I″_δh
end

function _step9(a, e, δe, e″δL, L, I, δI, sin_half_I″_δh, h, Σlgh)
    # Compute the osculating elements
    I1 = a
    I2 = e * cos(atan(e * sin(L), 1 + e * cos(L))) + δe
    I3 = e * sin(atan(e * sin(L), 1 + e * cos(L)))
    I4 = sin(I / 2 + δI / 2) * cos(h)
    I5 = sin(I / 2 + δI / 2) * sin(h)
    I6 = Σlgh

    return I1, I2, I3, I4, I5, I6
end

function _step10(a, h, k, p, q, L, aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, fⱼ₂, Lⱼ₂, γ, γ′, θ, rⱼ₂, η, νⱼ₂)
    # Inverse of step 5
    δa = -2 * γ * rⱼ₂ * (1 + (1 + η) * sin(νⱼ₂ + gⱼ₂)^2)
    a_osc = aⱼ₂ + δa

    # Inverse of step 6-9 to get the osculating ModEqN elements
    δh = γ′ * sin(νⱼ₂ + gⱼ₂) * cos(νⱼ₂ + gⱼ₂) / tan(Iⱼ₂)

    # Compute corrections
    δe = (2 * γ / (aⱼ₂ * η)) * (1 + η) * sin(νⱼ₂ + gⱼ₂) * cos(νⱼ₂ + gⱼ₂)
    e_osc = eⱼ₂ + δe

    δI = -2 * γ′ * cos(νⱼ₂ + gⱼ₂) * sin(νⱼ₂ + gⱼ₂) / sin(Iⱼ₂)
    I_osc = Iⱼ₂ + δI

    # Compute the longitude correction
    e″δL = -θ * ((2 + eⱼ₂^2) / eⱼ₂) * cos(2 * νⱼ₂ + 2 * gⱼ₂) - γ′ * cos(νⱼ₂ + gⱼ₂) / (eⱼ₂ * tan(Iⱼ₂))
    Σlgh = Lⱼ₂ + θ * (1 + 1 / (1 - eⱼ₂^2)) * sin(νⱼ₂ + gⱼ₂)^2 - δh
    L_osc = Σlgh - gⱼ₂ - hⱼ₂

    # Convert to ModEqN
    n_osc = sqrt(μ / a_osc^3)
    h_osc = e_osc * sin(gⱼ₂)
    k_osc = e_osc * cos(gⱼ₂)
    p_osc = tan(I_osc / 2) * sin(hⱼ₂)
    q_osc = tan(I_osc / 2) * cos(hⱼ₂)

    return n_osc, h_osc, k_osc, p_osc, q_osc, L_osc
end
