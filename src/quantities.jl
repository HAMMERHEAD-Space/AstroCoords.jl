export meanMotion
"""    meanMotion(a::Number, μ::Number)

Computes the Keplerian mean motion about a central body.

# Arguments
-`a::Number`: The semi-major axis of the orbit.
-`μ::Number`: Standard graviational parameter of central body.

# Returns
- `n::Number`: The orbital mean motion.
"""
function meanMotion(a::Number, μ::Number)
    return √(μ / (a^3.0))
end

"""    meanMotion(X::AstroCoord, μ::Number)

Computes the Keplerian mean motion about a central body.

# Arguments
-`X::AstroCoord`: An coordinate set describing the orbit.
-`μ::Number`: Standard graviational parameter of central body.

# Returns
- `n::Number`: The orbital mean motion.
"""
function meanMotion(X::AstroCoord, μ::Number, args...)
    kep = Keplerian(X, μ, args...)

    return meanMotion(kep.a, μ)
end

"""    meanMotion(X::EDromo, μ::Number)

Computes the Keplerian mean motion about a central body from EDromo coordinates.
This is an optimized version that extracts the semi-major axis directly from ζ₃.

# Arguments
-`X::EDromo`: An EDromo coordinate set describing the orbit.
-`μ::Number`: Standard graviational parameter of central body.

# Returns
- `n::Number`: The orbital mean motion.
"""
function meanMotion(X::EDromo, μ::Number, args...)
    # ζ₃ = -μ/(2E) = a (semi-major axis)
    return meanMotion(X.ζ₃, μ)
end

"""    meanMotion(X::StiefelScheifele, μ::Number)

Computes the Keplerian mean motion about a central body from Stiefel-Scheifele coordinates.
This is an optimized version that computes the semi-major axis directly from ω.

# Arguments
-`X::StiefelScheifele`: A Stiefel-Scheifele coordinate set describing the orbit.
-`μ::Number`: Standard graviational parameter of central body.

# Returns
- `n::Number`: The orbital mean motion.
"""
function meanMotion(X::StiefelScheifele, μ::Number, args...)
    # ω² = -E/2, and E = -μ/(2a), so a = μ/(2ω²)
    a = μ / (2.0 * X.ω^2)
    return meanMotion(a, μ)
end

export orbitalPeriod
"""    orbitalPeriod(a::Number, μ::Number)

Computes the Keplerian orbital period about a central body.

# Arguments
-`a::Number`: The semi-major axis of the orbit.
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`T::Number`: The orbital period.
"""
function orbitalPeriod(a::Number, μ::Number)
    return 2.0 * π / √(μ / (a^3.0))
end

"""    orbitalPeriod(X::AstroCoord, μ::Number)

Computes the Keplerian orbital period about a central body.

# Arguments
-`X::AstroCoord`: An coordinate set describing the orbit.
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`T::Number`: The orbital period.
"""
function orbitalPeriod(X::AstroCoord, μ::Number, args...)
    kep = Keplerian(X, μ, args...)

    return orbitalPeriod(kep.a, μ)
end

"""    orbitalPeriod(X::EDromo, μ::Number)

Computes the Keplerian orbital period about a central body from EDromo coordinates.
This is an optimized version that extracts the semi-major axis directly from ζ₃.

# Arguments
-`X::EDromo`: An EDromo coordinate set describing the orbit.
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`T::Number`: The orbital period.
"""
function orbitalPeriod(X::EDromo, μ::Number, args...)
    # ζ₃ = -μ/(2E) = a (semi-major axis)
    return orbitalPeriod(X.ζ₃, μ)
end

"""    orbitalPeriod(X::StiefelScheifele, μ::Number)

Computes the Keplerian orbital period about a central body from Stiefel-Scheifele coordinates.
This is an optimized version that computes the semi-major axis directly from ω.

# Arguments
-`X::StiefelScheifele`: A Stiefel-Scheifele coordinate set describing the orbit.
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`T::Number`: The orbital period.
"""
function orbitalPeriod(X::StiefelScheifele, μ::Number, args...)
    # ω² = -E/2, and E = -μ/(2a), so a = μ/(2ω²)
    a = μ / (2.0 * X.ω^2)
    return orbitalPeriod(a, μ)
end

export orbitalNRG
"""    orbitalNRG(a::Number, μ::Number)

Computes the keplerian orbital energy.

# Arguments
-`a::Number`: The semi-major axis of the orbit.
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`NRG::Number`: The orbital energy. 
"""
function orbitalNRG(a::Number, μ::Number)
    return -μ / (2.0 * a)
end

"""    orbitalNRG(X::AstroCoord, μ::Number)

Computes the keplerian orbital energy.

# Arguments
-`X::AstroCoord`: An coordinate set describing the orbit.
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`NRG::Number`: The orbital energy. 
"""
function orbitalNRG(X::AstroCoord, μ::Number, args...)
    kep = Keplerian(X, μ, args...)

    return orbitalNRG(kep.a, μ)
end

export angularMomentumVector
"""    angularMomentumVector(u::AbstractVector{<:Number})

Computes the instantaneous angular momentum vector from a Cartesian state vector.

# Arguments
-`u::AbstractVector{<:Number}`: The Cartesian state vector [x; y; z; ẋ; ẏ; ż].

# Returns
-'angular_momentum::Vector{<:Number}': 3-Dimensional angular momentum vector.
"""
function angularMomentumVector(u::AbstractVector{<:Number})
    r = SVector{3}(u[1], u[2], u[3])
    v = SVector{3}(u[4], u[5], u[6])

    return cross(r, v)
end

"""    angularMomentumVector(X::AstroCoord, μ::Number)

Computes the instantaneous angular momentum vector from a Cartesian state vector.

# Arguments
-`X::AstroCoord`: An coordinate set describing the orbit.
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`angular_momentum::Vector{<:Number}`: 3-Dimensional angular momentum vector.
"""
function angularMomentumVector(X::AstroCoord, μ::Number, args...)
    cart = Cartesian(X, μ, args...)

    return angularMomentumVector(params(cart))
end

export angularMomentumQuantity
"""    angularMomentumQuantity(u::AbstractVector{<:Number})

Computes the instantaneous angular momentum.

# Arguments
-`u::AbstractVector{<:Number}`: The Cartesian state vector [x; y; z; ẋ; ẏ; ż].

# Returns
-`angular_momentum::Number`: Angular momentum of the body.
"""
function angularMomentumQuantity(u::AbstractVector{<:Number})
    return norm(angularMomentumVector(u))
end

"""    angularMomentumQuantity(X::AstroCoord, μ::Number)

Computes the instantaneous angular momentum.

# Arguments
-`X::AstroCoord`: An coordinate set describing the orbit.
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`angular_momentum::Number`: Angular momentum of the body.
"""
function angularMomentumQuantity(X::AstroCoord, μ::Number, args...)
    cart = Cartesian(X, μ, args...)

    return angularMomentumQuantity(params(cart))
end
