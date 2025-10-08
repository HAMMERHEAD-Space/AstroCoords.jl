#TODO: ASTRO PROBLEMS/SOLVERS PKG??
#TODO: DEFINE ADJOINT (IMPLICIT?)
export KeplerSolver
"""
    KeplerSolver(M::T, e::Number; tol::Float64=10 * eps(T)) where {T<:Number}

Solves for true anomaly given the mean anomaly and eccentricity of an orbit.

# Arguments
-`M::Number`: Mean Anomaly of the orbit [radians].
-`e::Number`: Eccentricity of the orbit.

# Keyword Arguments
-`tol::Float64`: Convergence tolerance of Kepler solver. [Default=10*eps(T)]

# Returns
-`f::Number``: True Anomaly of the orbit [radians]
"""
function KeplerSolver(M::T, e::Number; tol::Float64=10 * eps(T)) where {T<:Number}
    if e < 1.0
        E_guess = M
        fE = M - (E_guess - e * sin(E_guess))
        dfE = (e * cos(E_guess) - 1)

        while abs(fE) > tol
            E_guess -= fE / dfE
            fE = M - (E_guess - e * sin(E_guess))
            dfE = (e * cos(E_guess) - 1)
        end
        E_guess = rem2pi(E_guess, RoundDown)

        # Compute square roots once to avoid repeated allocations
        E_half = E_guess / 2
        sin_E_half = sin(E_half)
        cos_E_half = cos(E_half)
        sqrt_1_plus_e = √(1 + e)
        sqrt_1_minus_e = √(1 - e)
        
        f = 2 * atan(sqrt_1_plus_e * sin_E_half, sqrt_1_minus_e * cos_E_half)

        if f < 0.0
            f = 2.0 * π - abs(f)
        end
    else
        F_guess = M
        fF = M + (F_guess - e * sinh(F_guess))
        dfF = (1.0 - e * cosh(F_guess))

        while abs(fF) > tol
            F_guess -= fF / dfF
            fF = M + (F_guess - e * sinh(F_guess))
            dfF = (1.0 - e * cosh(F_guess))
        end

        F_guess = rem2pi(F_guess, RoundDown)

        # Compute square roots once to avoid repeated allocations
        F_half = F_guess / 2
        sinh_F_half = sinh(F_half)
        cosh_F_half = cosh(F_half)
        sqrt_1_plus_e = √(1 + e)
        sqrt_e_minus_1 = √(e - 1)
        
        f = 2.0 * atan(sqrt_1_plus_e * sinh_F_half, sqrt_e_minus_1 * cosh_F_half)

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
        F = atanh(sin(f) * √(e^2 - 1.0) / (1.0 + e * cos(f)))

        F = rem2pi(F, RoundDown)

        E = F
    end

    return E
end

export eccentricAnomaly2MeanAnomaly
"""
    eccentricAnomaly2MeanAnomaly(E::Number, e::Number)

Converts the eccentric anomaly into the mean anomaly.

# Arguments
-`E::Number`: Eccentric anomaly of the orbit [radians].
-`e::Number`: Eccentricity of the orbit.

# Returns
-`M::Number`: Mean anomaly of the orbit [radians].
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

export meanAnomaly2TrueAnomaly
"""
    meanAnomaly2TrueAnomaly(M::Number, e::Number)

Converts the mean anomaly into the true anomaly.

# Arguments
-`M::Number`: Mean anomaly of the orbit [radians].
-`e::Number`: Eccentricity of the orbit.

# Returns
-`f::Number`: True anomaly of the orbit [radians].
"""
@inline function meanAnomaly2TrueAnomaly(M::Number, e::Number)
    f = KeplerSolver(M, e)

    return f
end
