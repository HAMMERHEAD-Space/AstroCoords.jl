#TODO: ASTRO PROBLEMS/SOLVERS PKG??
#TODO: DEFINE ADJOINT (IMPLICIT?)
export KeplerSolver
"""
    KeplerSolver(M::T, e::Number; tol::Float64=10 * eps(T), max_iters::Int=100) where {T<:Number}

Solves for true anomaly given the mean anomaly and eccentricity of an orbit.

# Arguments
-`M::Number`: Mean Anomaly of the orbit [radians].
-`e::Number`: Eccentricity of the orbit.

# Keyword Arguments
-`tol::Float64`: Convergence tolerance of Kepler solver. [Default=10*eps(T)]
-`max_iters::Int`: Maximum number of iterations. [Default=100]

# Returns
-`f::Number``: True Anomaly of the orbit [radians]
"""
function KeplerSolver(
    M::T, e::Number; tol::Float64=10 * eps(T), max_iters::Int=100
) where {T<:Number}
    if e < 1.0
        E_guess = M
        fE = M - (E_guess - e * sin(E_guess))
        dfE = (e * cos(E_guess) - 1)

        iter = 0
        while abs(fE) > tol && iter < max_iters
            E_guess -= fE / dfE
            fE = M - (E_guess - e * sin(E_guess))
            dfE = (e * cos(E_guess) - 1)
            iter += 1
        end
        E_guess = rem2pi(E_guess, RoundDown)

        f = 2 * atan(√(1 + e) * sin(E_guess / 2), √(1 - e) * cos(E_guess / 2))

        if f < 0.0
            f = 2.0 * π - abs(f)
        end
    else
        F_guess = M
        fF = M + (F_guess - e * sinh(F_guess))
        dfF = (1.0 - e * cosh(F_guess))

        iter = 0
        while abs(fF) > tol && iter < max_iters
            F_guess -= fF / dfF
            fF = M + (F_guess - e * sinh(F_guess))
            dfF = (1.0 - e * cosh(F_guess))
            iter += 1
        end

        F_guess = rem2pi(F_guess, RoundDown)

        f = 2.0 * atan(√(1 + e) * sinh(F_guess / 2), √(e - 1) * cosh(F_guess / 2))

        if f < 0.0
            f = 2.0 * π - abs(f)
        end
    end

    return f
end

export trueAnomaly2MeanAnomaly
"""
    trueAnomaly2MeanAnomaly(f::Number, e::Number)

Converts the true anomaly into the mean anomaly.

# Arguments
-`f::Number`: True anomaly of the orbit [radians].
-`e::Number`: Eccentricity of the orbit.

# Returns
-`M::Number`: Mean anomaly of the orbit [radians].
"""
@inline function trueAnomaly2MeanAnomaly(f::Number, e::Number)
    E = trueAnomaly2EccentricAnomaly(f, e)
    M = eccentricAnomaly2MeanAnomaly(E, e)

    return M
end

export trueAnomaly2EccentricAnomaly
"""
    trueAnomaly2EccentricAnomaly(f::Number, e::Number)

Converts the true anomaly into the mean anomaly.

# Arguments
-`f::Number`: True anomaly of the orbit [radians].
-`e::Number`: Eccentricity of the orbit.

# Returns
-`E::Number`: Eccentric anomaly of the orbit [radians].
"""
@inline function trueAnomaly2EccentricAnomaly(f::Number, e::Number)
    if e < 1.0
        E = atan(
            (sin(f) * √(1 - e^2)) / (1.0 + e * cos(f)), (e + cos(f)) / (1.0 + e * cos(f))
        )

        E = rem2pi(E, RoundDown)
    else
        E = 2.0 * atanh(√((e - 1.0) / (1.0 + e)) * tan(f / 2.0))

        E = rem2pi(E, RoundDown)
    end

    return E
end

export eccentricAnomaly2MeanAnomaly
"""
    eccentricAnomaly2MeanAnomaly(E::Number, e::Number)

Converts the true anomaly into the mean anomaly.

# Arguments
-`E::Number`: Eccentric anomaly of the orbit [radians].
-`e::Number`: Eccentricity of the orbit.

# Returns
-'M::Number': Mean anomaly of the orbit [radians].
"""
@inline function eccentricAnomaly2MeanAnomaly(E::Number, e::Number)
    if e < 1.0
        M = E - e * sin(E)

        M = rem2pi(M, RoundDown)
    else
        M = e * sinh(E) - E

        M = rem2pi(M, RoundDown)
    end

    return M
end

export eccentricAnomaly2TrueAnomaly
"""
    eccentricAnomaly2TrueAnomaly(E::Number, e::Number)

Converts the eccentric anomaly into the true anomaly.

# Arguments
-`E::Number`: Eccentric anomaly of the orbit [radians].
-`e::Number`: Eccentricity of the orbit.

# Returns
-`f::Number`: True anomaly of the orbit [radians].
"""
@inline function eccentricAnomaly2TrueAnomaly(E::Number, e::Number)
    if e < 1.0
        f = 2.0 * atan(√(1.0 + e) * sin(E / 2.0), √(1.0 - e) * cos(E / 2.0))
    else
        f = atan(√(e + 1.0) * sinh(E / 2.0), √(e - 1.0) * cosh(E / 2.0))
    end
    return f
end

export meanAnomaly2TrueAnomaly
"""
    meanAnomaly2TrueAnomaly(M::T, e::Number; tol::Float64=10 * eps(T), max_iters::Int=100) where {T<:Number}

Converts the mean anomaly into the true anomaly.

# Arguments
-`M::Number`: Mean anomaly of the orbit [radians].
-`e::Number`: Eccentricity of the orbit.

# Keyword Arguments
-`tol::Float64`: Convergence tolerance of Kepler solver. [Default=10*eps(T)]
-`max_iters::Int`: Maximum number of iterations. [Default=100]

# Returns
-`f::Number`: Mean anomaly of the orbit [radians].
"""
@inline function meanAnomaly2TrueAnomaly(
    M::T, e::Number; tol::Float64=10 * eps(T), max_iters::Int=100
) where {T<:Number}
    return KeplerSolver(M, e; tol=tol, max_iters=max_iters)
end

export meanAnomaly2EccentricAnomaly
"""
    meanAnomaly2EccentricAnomaly(M::T, e::Number; tol::Float64=10 * eps(T), max_iters::Int=100) where {T<:Number}

Converts the Mean Anomaly into the Eccentric Anomaly

# Arguments
-`M::Number`: Mean Anomaly of the orbit [radians]
-`e::Number`: Eccentricity of the orbit

# Keyword Arguments
-`tol::Float64`: Convergence tolerance of Kepler solver. [Default=10*eps(T)]
-`max_iters::Int`: Maximum number of iterations. [Default=100]

# Returns
-`E::Number`: Eccentric Anomaly of the orbit [radians]
"""
@inline function meanAnomaly2EccentricAnomaly(
    M::T, e::Number; tol::Float64=10 * eps(T), max_iters::Int=100
) where {T<:Number}
    return trueAnomaly2EccentricAnomaly(
        meanAnomaly2TrueAnomaly(M, e; tol=tol, max_iters=max_iters), e
    )
end
