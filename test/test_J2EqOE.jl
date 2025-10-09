@testset "J2EqOE Conversion" begin
    μ = 3.986004415e+05
    R = 6.378137e+03
    J2 = 1.0826261738522e-03

    n″ = 1.0487861550167e-03
    h″ = 1.4260897234681e-03
    k″ = -9.2912382513865e-03
    p″ = 6.6371635585001e-01
    q″ = -3.2380704393120e-01
    L″ = 4.8735182761240e+00

    u = [n″, h″, k″, p″, q″, L″]

    u_J2IOE = AstroCoords.modEqN2IOE(u, μ)

    expected_J2IOE = [
        7.1294270251742e+03,
        8.9420945754925e-03,
        -2.8982382142187e-03,
        5.3390774990364e-01,
        -2.6047736914253e-01,
        4.8735182761240e+00,
    ]

    @test u_J2IOE ≈ expected_J2IOE rtol = 1e-13

    @testset "Step 1" begin
        J2IOE = [
            7.1294270251742e+03,
            8.9420945754925e-03,
            -2.8982382142187e-03,
            5.3390774990364e-01,
            -2.6047736914253e-01,
            4.8735182761240e+00,
        ]
        aⱼ₂, eⱼ₂, Iⱼ₂, hⱼ₂, gⱼ₂, f, Lⱼ₂ = AstroCoords._step1(J2IOE, μ)

        expected_aⱼ₂ = 7.1294270251742e+03
        expected_eⱼ₂ = 9.4000446883730e-03
        expected_Iⱼ₂ = 1.2721904014791e+00
        expected_hⱼ₂ = 2.0246926216394e+00
        expected_gⱼ₂ = 9.6460100038126e-01
        expected_Lⱼ₂ = 1.8842246541034e+00

        @test aⱼ₂ ≈ expected_aⱼ₂ rtol = 1e-13
        @test eⱼ₂ ≈ expected_eⱼ₂ rtol = 1e-13
        @test hⱼ₂ ≈ expected_hⱼ₂ rtol = 1e-13
        @test Iⱼ₂ ≈ expected_Iⱼ₂ rtol = 1e-13
        @test gⱼ₂ ≈ expected_gⱼ₂ rtol = 1e-13
        @test Lⱼ₂ ≈ expected_Lⱼ₂ rtol = 1e-13
    end

    @testset "Step 2" begin
        aⱼ₂ = 7.1294270251742e+03
        eⱼ₂ = 9.4000446883730e-03
        Iⱼ₂ = 1.2721904014791e+00

        k, η, γ, γ′, θ = AstroCoords._step2(J2, R, eⱼ₂, aⱼ₂, Iⱼ₂)

        expected_k = 2.2020958264503e+04
        expected_η = 9.9995581860393e-01
        expected_γ = 4.3323841440301e-04
        expected_γ′ = 4.3331498717247e-04
        expected_θ = 2.9418810951484e-01

        @test k ≈ expected_k rtol = 1e-13
        @test η ≈ expected_η rtol = 1e-13
        @test γ ≈ expected_γ rtol = 1e-13
        @test γ′ ≈ expected_γ′ rtol = 1e-13
        @test θ ≈ expected_θ rtol = 1e-13
    end

    @testset "Step 3" begin
        fⱼ₂ = 1.9020433345198198
        eⱼ₂ = 9.4000446883730e-03

        Eⱼ₂ = AstroCoords._step3(fⱼ₂, eⱼ₂)

        expected_Eⱼ₂ = 1.8931405532060e+00

        @test Eⱼ₂ ≈ expected_Eⱼ₂ rtol = 1e-13
    end

    @testset "Step 4" begin
        aⱼ₂ = 7.1294270251742e+03
        eⱼ₂ = 9.4000446883730e-03
        η = 9.9995581860393e-01
        Eⱼ₂ = 1.8931405532060e+00
        Lⱼ₂ = 1.8842246541034e+00

        rⱼ₂, νⱼ₂, _ = AstroCoords._step4(aⱼ₂, eⱼ₂, η, Eⱼ₂, Lⱼ₂)

        expected_rⱼ₂ = 7.1506573807183e+03
        expected_νⱼ₂ = 1.9020433345198e+00

        @test rⱼ₂ ≈ expected_rⱼ₂ rtol = 1e-13
        @test νⱼ₂ ≈ expected_νⱼ₂ rtol = 1e-13
    end

    @testset "Step 5" begin
        aⱼ₂ = 7.1294270251742e+03
        γ = 4.3323841440301e-04
        γ′ = 4.3331498717247e-04
        θ = 2.9418810951484e-01
        rⱼ₂ = 7.1506573807183e+03
        η = 9.9995581860393e-01
        eⱼ₂ = 9.4000446883730e-03
        gⱼ₂ = 9.6460100038126e-01
        νⱼ₂ = 1.9020433345198e+00
        Lⱼ₂ = 1.8842246541034e+00

        a, δh = AstroCoords._step5(aⱼ₂, γ, γ′, θ, rⱼ₂, η, eⱼ₂, gⱼ₂, νⱼ₂, Lⱼ₂)

        expected_a = 7.1366000000000e+03
        expected_δh = -1.1070091328349e-04

        @test a ≈ expected_a rtol = 1e-13
        @test δh ≈ expected_δh rtol = 1e-13
    end

    @testset "Step 6" begin
        γ′ = 4.3331498717247e-04
        θ = 2.9418810951484e-01
        νⱼ₂ = 1.9020433345198e+00
        Lⱼ₂ = 1.8842246541034e+00
        eⱼ₂ = 9.4000446883730e-03
        gⱼ₂ = 9.6460100038126e-01
        hⱼ₂ = 2.0246926216394e+00
        δh = -1.1070091328349e-04

        Σlgh = AstroCoords._step6(γ′, θ, νⱼ₂, Lⱼ₂, eⱼ₂, gⱼ₂, hⱼ₂, δh)

        expected_Σlgh = 4.8729592715682e+00

        @test Σlgh ≈ expected_Σlgh rtol = 1e-13
    end

    @testset "Step 7" begin
        νⱼ₂ = 1.9020433345198e+00
        eⱼ₂ = 9.4000446883730e-03
        η = 9.9995581860393e-01
        aⱼ₂ = 7.1294270251742e+03
        rⱼ₂ = 7.1506573807183e+03
        γ = 4.3323841440301e-04
        γ′ = 4.3331498717247e-04
        θ = 2.9418810951484e-01
        gⱼ₂ = 9.6460100038126e-01

        v1, v2, v3, v4, δe, e″δL = AstroCoords._step7(νⱼ₂, eⱼ₂, η, aⱼ₂, rⱼ₂, γ, γ′, θ, gⱼ₂)

        expected_v1 = -9.7268781528180e-01
        expected_v2 = -9.5884220957866e-01
        expected_v3 = -9.6354316647622e-01
        expected_v4 = 1.9910139555114e+00
        expected_δe = 8.1222972039185e-05
        expected_e″δL = -4.0701787616281e-04

        @test v1 ≈ expected_v1 rtol = 1e-13
        @test v2 ≈ expected_v2 rtol = 1e-13
        @test v3 ≈ expected_v3 rtol = 1e-13
        @test v4 ≈ expected_v4 rtol = 1e-13
        @test δe ≈ expected_δe rtol = 1e-12
        @test e″δL ≈ expected_e″δL rtol = 1e-13
    end

    @testset "Step 8" begin
        γ′ = 4.3331498717247e-04
        θ = 2.9418810951484e-01
        νⱼ₂ = 1.9020433345198e+00
        eⱼ₂ = 9.4000446883730e-03
        gⱼ₂ = 9.6460100038126e-01
        δh = -1.1070091328349e-04
        Iⱼ₂ = 1.2721904014791e+00

        δI, sin_half_I″_δh = AstroCoords._step8(γ′, θ, νⱼ₂, eⱼ₂, gⱼ₂, δh, Iⱼ₂)

        expected_δI = 1.5460976151790e-04
        expected_sin_half_I″_δh = -6.5762859846059e-05

        @test δI ≈ expected_δI rtol = 1e-13
        @test sin_half_I″_δh ≈ expected_sin_half_I″_δh rtol = 1e-13
    end

    @testset "Step 9" begin
        a = 7.1366000000000e+03
        eⱼ₂ = 9.4000446883730e-03
        δe = 8.1222972039185e-05
        e″δL = -4.0701787616281e-04
        Lⱼ₂ = 1.8842246541034e+00
        Iⱼ₂ = 1.2721904014791e+00
        δI = 1.5460976151790e-04
        sin_half_I″_δh = -6.5762859846059e-05
        hⱼ₂ = 2.0246926216394e+00
        Σlgh = 4.8729592715682e+00

        I1, I2, I3, I4, I5, I6 = AstroCoords._step9(
            a, eⱼ₂, δe, e″δL, Lⱼ₂, Iⱼ₂, δI, sin_half_I″_δh, hⱼ₂, Σlgh
        )

        expected_I1 = 7.1366000000000e+03
        expected_I2 = 9.4812669855769e-03
        expected_I3 = 1.2737365091406e+00
        expected_I4 = 2.0245275207411e+00
        expected_I5 = 9.6460100038126e-01
        expected_I6 = 4.8729592715682e+00

        @test I1 ≈ expected_I1 rtol = 1e-13
        @test I2 ≈ expected_I2 rtol = 1e-13
        @test I3 ≈ expected_I3 rtol = 1e-13
        @test I4 ≈ expected_I4 rtol = 1e-13
        @test I5 ≈ expected_I5 rtol = 1e-13
        @test I6 ≈ expected_I6 rtol = 1e-13
    end
end
