export meanMotion
"""
    meanMotion(a::Number, μ::Number)

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

"""
    meanMotion(X::AstroCoord, μ::Number)

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

export orbitalPeriod
"""
    orbitalPeriod(a::Number, μ::Number)

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

"""
    orbitalPeriod(X::AstroCoord, μ::Number)

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

export orbitalNRG
"""
    orbitalNRG(a::Number, μ::Number)

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

"""
    orbitalNRG(X::AstroCoord, μ::Number)

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
"""
    angularMomentumVector(u::AbstractVector{<:Number})

Computes the instantaneous angular momentum vector from a Cartesian state vector.

# Arguments
-`u::AbstractVector{<:Number}`: The Cartesian state vector [x; y; z; ẋ; ẏ; ż].

# Returns
-'angular_momentum::Vector{<:Number}': 3-Dimensional angular momemtum vector.
"""
function angularMomentumVector(u::AbstractVector{<:Number})
    r = SVector{3}(u[1], u[2], u[3])
    v = SVector{3}(u[4], u[5], u[6])

    return cross(r, v)
end

"""
    angularMomentumVector(X::AstroCoord, μ::Number)

Computes the instantaneous angular momentum vector from a Cartesian state vector.

# Arguments
-`X::AstroCoord`: An coordinate set describing the orbit.
-`μ::Number`: Standard graviational parameter of central body.

# Returns
-`angular_momentum::Vector{<:Number}`: 3-Dimensional angular momemtum vector.
"""
function angularMomentumVector(X::AstroCoord, μ::Number, args...)
    cart = Cartesian(X, μ, args...)

    return angularMomentumVector(params(cart))
end

export angularMomentumQuantity
"""
    angularMomentumQuantity(u::AbstractVector{<:Number})

Computes the instantaneous angular momentum.

# Arguments
-`u::AbstractVector{<:Number}`: The Cartesian state vector [x; y; z; ẋ; ẏ; ż].

# Returns
-`angular_momentum::Number`: Angular momentum of the body.
"""
function angularMomentumQuantity(u::AbstractVector{<:Number})
    return norm(angularMomentumVector(u))
end

"""
    angularMomentumQuantity(X::AstroCoord, μ::Number)

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
