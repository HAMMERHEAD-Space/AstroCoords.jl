@testset "Test Anomaly Conversions" begin
    M = 2π * rand()
    e = rand()

    E = meanAnomaly2EccentricAnomaly(M, e)
    f = meanAnomaly2TrueAnomaly(M, e)

    E2 = trueAnomaly2EccentricAnomaly(f, e)
    f2 = eccentricAnomaly2TrueAnomaly(E, e)

    @test E ≈ E2 atol = 1e-14
    @test f ≈ f2 atol = 1e-14

    M_from_E = eccentricAnomaly2MeanAnomaly(E, e)
    M_from_f = trueAnomaly2MeanAnomaly(f, e)

    @test M_from_E ≈ M atol = 1e-14
    @test M_from_f ≈ M atol = 1e-14
end
