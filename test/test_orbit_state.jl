@testset "OrbitState" begin
    u0 = [-1076.2, -6765.9, -332.3, 9.357, -3.312, -1.188]
    μ = 398600.4415
    ep = 0.0
    frm = :ICRF

    s = OrbitState(Cartesian(u0), ep, frm)

    @testset "construction and accessors" begin
        @test coords(s) == Cartesian(u0)
        @test epoch(s) == ep
        @test frame(s) == frm
    end

    @testset "construction type" begin
        @test s isa OrbitState{<:Cartesian,<:Any,:ICRF}
        @test frame(s) === :ICRF
    end

    @testset "convert_coords Cartesian → Keplerian" begin
        s_kep = convert_coords(s, Keplerian, μ)
        @test s_kep isa OrbitState{<:Keplerian,<:Any,:ICRF}
        @test epoch(s_kep) == ep
        @test frame(s_kep) == frm
    end

    @testset "convert_coords round-trip" begin
        s_kep = convert_coords(s, Keplerian, μ)
        s_rt = convert_coords(s_kep, Cartesian, μ)
        @test s_rt isa OrbitState{<:Cartesian,<:Any,:ICRF}
        @test epoch(s_rt) == ep
        @test frame(s_rt) == frm
        @test params(coords(s_rt)) ≈ params(coords(s)) rtol = 1e-10
    end

    @testset "change_frame Cartesian" begin
        frames = FrameSystem{2,Float64}()
        add_axes!(frames, :ICRF, 1)
        add_axes_fixed_angles!(frames, :ITRF, 2, 1, [0.0, 0.0, 0.0], :ZYX)
        cr_fwd = compile_rotation6(frames, :ICRF, :ITRF)
        cr_inv = compile_rotation6(frames, :ITRF, :ICRF)

        s_itrf = change_frame(s, Val{:ITRF}(), cr_fwd, μ)
        @test s_itrf isa OrbitState{<:Cartesian,<:Any,:ITRF}
        @test epoch(s_itrf) == ep
        @test frame(s_itrf) === :ITRF
        # identity rotation → coords unchanged
        @test params(coords(s_itrf)) ≈ params(coords(s))

        # round-trip
        s_rt = change_frame(s_itrf, Val{:ICRF}(), cr_inv, μ)
        @test s_rt isa OrbitState{<:Cartesian,<:Any,:ICRF}
        @test frame(s_rt) === :ICRF
        @test params(coords(s_rt)) ≈ params(coords(s))
    end

    @testset "change_frame Symbol convenience" begin
        frames = FrameSystem{2,Float64}()
        add_axes!(frames, :ICRF, 1)
        add_axes_fixed_angles!(frames, :ITRF, 2, 1, [0.0, 0.0, 0.0], :ZYX)

        s_itrf = change_frame(s, :ITRF, frames, μ)
        @test s_itrf isa OrbitState{<:Cartesian,<:Any,:ITRF}
        @test frame(s_itrf) === :ITRF
    end

    @testset "change_frame Keplerian" begin
        frames = FrameSystem{2,Float64}()
        add_axes!(frames, :ICRF, 1)
        add_axes_fixed_angles!(frames, :ITRF, 2, 1, [0.0, 0.0, 0.0], :ZYX)
        cr_fwd = compile_rotation6(frames, :ICRF, :ITRF)

        s_kep = convert_coords(s, Keplerian, μ)
        s_kep_itrf = change_frame(s_kep, Val{:ITRF}(), cr_fwd, μ)
        @test s_kep_itrf isa OrbitState{<:Keplerian,<:Any,:ITRF}
        @test epoch(s_kep_itrf) == ep
        @test frame(s_kep_itrf) === :ITRF
        # identity rotation → Keplerian elements preserved
        @test params(coords(s_kep_itrf)) ≈ params(coords(s_kep))
    end

    @testset "Tempo.Epoch support" begin
        using Tempo
        ep_tempo = Epoch(0.0, TDB)
        s_tempo = OrbitState(Cartesian(u0), ep_tempo, frm)
        @test epoch(s_tempo) === ep_tempo
        s_kep = convert_coords(s_tempo, Keplerian, μ)
        @test epoch(s_kep) === ep_tempo
        @test frame(s_kep) == frm
    end
end
