########################################################################################
#
# Tests for Differentiation Across the Coordinate Sets
#
########################################################################################
# Currently Supported & Tested
# Diffractor, Enzyme, ForwardDiff, FiniteDiff, Mooncake, PolyesterForwardDiff, Zygote
########################################################################################
const _BACKENDS = (
    ("ForwardDiff", AutoForwardDiff()),
    ("Diffractor", AutoDiffractor()),
    ("Enzyme", AutoEnzyme()),
    ("Mooncake", AutoMooncake(; config=nothing)),
    ("PolyesterForwardDiff", AutoPolyesterForwardDiff()),
)

@testset "Coordinate Transformation Differentiation" begin
    state = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ]

    μ = 3.986004415e5

    for backend in _BACKENDS
        testset_name = "Coordinate Set Transformation " * backend[1]
        @testset "$testset_name" begin
            for set in _COORDINATE_SETS
                f_fd, df_fd = value_and_jacobian(
                    (x) -> set(Cartesian(x), μ), AutoFiniteDiff(), state
                )

                f_fd2, df_fd2 = value_and_derivative(
                    (x) -> set(Cartesian(state), x), AutoFiniteDiff(), μ
                )

                f_ad, df_ad = value_and_jacobian(
                    (x) -> Array(params(set(Cartesian(x), μ))), backend[2], state
                )

                @test f_fd == f_ad

                #TODO: Diffractor has some issue with the USM sets
                if backend[1] == "Diffractor"
                    @test df_fd ≈ df_ad rtol = 2e0
                else
                    @test df_fd ≈ df_ad atol = 1e-2
                end

                f_ad2, df_ad2 = value_and_derivative(
                    (x) -> Array(params(set(Cartesian(state), x))), backend[2], μ
                )

                @test f_fd2 == f_ad2
                @test df_fd2 ≈ df_ad2 atol = 1e-4

                coord = Array(params(set(Cartesian(state), μ)))

                f_fd3, df_fd3 = value_and_jacobian(
                    (x) -> Cartesian(set(x), μ), AutoFiniteDiff(), coord
                )

                f_fd4, df_fd4 = value_and_derivative(
                    (x) -> Cartesian(set(coord), x), AutoFiniteDiff(), μ
                )

                f_ad3, df_ad3 = value_and_jacobian(
                    (x) -> Cartesian(set(x), μ), backend[2], coord
                )

                f_ad4, df_ad4 = value_and_derivative(
                    (x) -> Cartesian(set(coord), x), backend[2], μ
                )

                @test f_fd3 == f_ad3

                if backend[1] == "Diffractor"
                    @test df_fd3 ≈ df_ad3 rtol = 2e0
                else
                    @test df_fd3 ≈ df_ad3 atol = 1e-2
                end

                @test f_fd4 == f_ad4
                @test df_fd4 ≈ df_ad4 atol = 1e-4
            end
        end
    end

    @testset "Coordinate Set Transformation Zygote" begin
        for set in _COORDINATE_SETS
            f_fd, df_fd = value_and_jacobian(
                (x) -> set(Cartesian(x), μ), AutoFiniteDiff(), state
            )

            f_fd2, df_fd2 = value_and_derivative(
                (x) -> set(Cartesian(state), x), AutoFiniteDiff(), μ
            )

            try
                f_ad, df_ad = value_and_jacobian(
                    (x) -> set(Cartesian(x), μ), AutoZygote(), state
                )
                @test f_fd == f_ad
                @test df_fd ≈ df_ad atol = 1e-2
            catch err
                @test err isa MethodError
                @test startswith(
                    sprint(showerror, err),
                    "MethodError: no method matching iterate(::Nothing)",
                )
            end

            try
                f_ad2, df_ad2 = value_and_derivative(
                    (x) -> set(Cartesian(state), x), AutoZygote(), μ
                )
                @test f_fd2 == f_ad2
                @test df_fd2 ≈ something.(df_ad2, 0.0) atol = 1e-4
            catch err
                @test err isa MethodError
                @test startswith(
                    sprint(showerror, err),
                    "MethodError: no method matching iterate(::Nothing)",
                )
            end

            coord = Array(params(set(Cartesian(state), μ)))

            f_fd3, df_fd3 = value_and_jacobian(
                (x) -> Cartesian(set(x), μ), AutoFiniteDiff(), coord
            )

            f_fd4, df_fd4 = value_and_derivative(
                (x) -> Cartesian(set(coord), x), AutoFiniteDiff(), μ
            )

            try
                f_ad3, df_ad3 = value_and_jacobian(
                    (x) -> Cartesian(set(x), μ), AutoZygote(), coord
                )
                @test f_fd3 == f_ad3
                @test df_fd3 ≈ df_ad3 atol = 1e-2
            catch err
                @test err isa MethodError
                @test startswith(
                    sprint(showerror, err),
                    "MethodError: no method matching iterate(::Nothing)",
                )
            end

            try
                f_ad4, df_ad4 = value_and_derivative(
                    (x) -> Cartesian(set(coord), x), AutoZygote(), μ
                )
                @test f_fd4 == f_ad4
                @test df_fd4 ≈ df_ad4 atol = 1e-4
            catch err
                @test err isa MethodError
                @test startswith(
                    sprint(showerror, err),
                    "MethodError: no method matching iterate(::Nothing)",
                )
            end
        end
    end
end
