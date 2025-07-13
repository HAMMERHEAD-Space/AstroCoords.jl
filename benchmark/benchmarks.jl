using AstroCoords
using BenchmarkTools

const SUITE = BenchmarkGroup()

SUITE["transformation"] = BenchmarkGroup(["convert"])
SUITE["anomalies"] = BenchmarkGroup(["keplerian math"])
SUITE["quantities"] = BenchmarkGroup(["derived"])

const _COORDINATE_SETS = [
    Cartesian,
    Delaunay,
    Keplerian,
    Milankovich,
    ModEq,
    Cylindrical,
    Spherical,
    USM7,
    USM6,
    USMEM,
    J2EqOE,
    EDromo,
    KustaanheimoStiefel,
]

const _state = [
    -1076.225324679696
    -6765.896364327722
    -332.3087833503755
    9.356857417032581
    -3.3123476319597557
    -1.1880157328553503
]

const _μ = 3.986004415e5

const _cart_state = Cartesian(_state)

for set in _COORDINATE_SETS
    SUITE["transformation"][string(set)] = @benchmarkable $(set)($_cart_state, $_μ)
    new_coord = set(_cart_state, _μ)
    SUITE["transformation"][string(set) * "reverse"] = @benchmarkable Cartesian(
        $new_coord, $_μ
    )
end

const _anomaly_conversions = [
    meanAnomaly2EccentricAnomaly,
    meanAnomaly2TrueAnomaly,
    trueAnomaly2EccentricAnomaly,
    trueAnomaly2MeanAnomaly,
    eccentricAnomaly2MeanAnomaly,
    eccentricAnomaly2TrueAnomaly,
]

const _e = 0.30230575359641376
const _angle = 0.6441434680007933

for f in _anomaly_conversions
    SUITE["anomalies"][string(f)] = @benchmarkable $(f)($_angle, $_e)
end

const _quantity_functions = [
    meanMotion, orbitalPeriod, orbitalNRG, angularMomentumVector, angularMomentumQuantity
]

for f in _quantity_functions
    SUITE["quantities"][string(f)] = @benchmarkable $(f)($_cart_state, $_μ)
end

# If a cache of tuned parameters already exists, use it, otherwise, tune and cache
# the benchmark parameters. Reusing cached parameters is faster and more reliable
# than re-tuning `suite` every time the file is included.
paramspath = joinpath(dirname(@__FILE__), "params.json")

if isfile(paramspath)
    loadparams!(SUITE, BenchmarkTools.load(paramspath)[1], :evals)
else
    tune!(SUITE)
    BenchmarkTools.save(paramspath, BenchmarkTools.params(SUITE))
end
