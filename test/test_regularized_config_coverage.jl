using Test
using AstroCoords
using StaticArrays
using Dates

@testset "RegularizedCoordinateConfig Coverage" begin
    @testset "Keyword Constructor" begin
        # Test default keyword constructor
        config = RegularizedCoordinateConfig()
        @test config.DU == 0.0
        @test config.TU == 0.0
        @test config.W == 0.0
        @test config.t₀ == 0.0
        @test config.flag_time isa PhysicalTime

        # Test with all keywords specified
        config = RegularizedCoordinateConfig(
            DU=1000.0,
            TU=60.0,
            W=0.5,
            t₀=100.0,
            flag_time=LinearTime()
        )
        @test config.DU == 1000.0
        @test config.TU == 60.0
        @test config.W == 0.5
        @test config.t₀ == 100.0
        @test config.flag_time isa LinearTime

        # Test with partial keywords
        config = RegularizedCoordinateConfig(DU=500.0, TU=30.0)
        @test config.DU == 500.0
        @test config.TU == 30.0
        @test config.W == 0.0
        @test config.t₀ == 0.0
        @test config.flag_time isa PhysicalTime
    end

    @testset "Positional Constructor" begin
        # Test direct positional constructor
        config = RegularizedCoordinateConfig(1500.0, 90.0, 0.3, 50.0, ConstantTime())
        @test config.DU == 1500.0
        @test config.TU == 90.0
        @test config.W == 0.3
        @test config.t₀ == 50.0
        @test config.flag_time isa ConstantTime
    end

    @testset "Type Promotion" begin
        # Test type promotion with mixed numeric types
        config = RegularizedCoordinateConfig(1000, 60.0, 0.5f0, 100, PhysicalTime())
        @test config.DU isa Number
        @test config.TU isa Number
        @test config.W isa Number
        @test config.t₀ isa Number
    end

    @testset "State-based Constructor" begin
        # Test constructor that computes DU and TU from state
        μ = 398600.4418  # Earth μ in km³/s²
        state = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]  # Circular orbit state

        config = RegularizedCoordinateConfig(state, μ)
        @test config.DU ≈ 7000.0 atol=1e-10
        @test config.TU ≈ sqrt(7000.0^3 / μ) atol=1e-6
        @test config.W == 0.0
        @test config.t₀ == 0.0
        @test config.flag_time isa PhysicalTime

        # Test with additional keyword arguments
        config = RegularizedCoordinateConfig(
            state, μ,
            W=0.2,
            t₀=200.0,
            flag_time=LinearTime()
        )
        @test config.DU ≈ 7000.0 atol=1e-10
        @test config.TU ≈ sqrt(7000.0^3 / μ) atol=1e-6
        @test config.W == 0.2
        @test config.t₀ == 200.0
        @test config.flag_time isa LinearTime
    end

    @testset "Flag Time Variants" begin
        # Test all three flag_time types
        config_physical = RegularizedCoordinateConfig(
            DU=1000.0, TU=60.0, flag_time=PhysicalTime()
        )
        @test config_physical.flag_time isa PhysicalTime

        config_constant = RegularizedCoordinateConfig(
            DU=1000.0, TU=60.0, flag_time=ConstantTime()
        )
        @test config_constant.flag_time isa ConstantTime

        config_linear = RegularizedCoordinateConfig(
            DU=1000.0, TU=60.0, flag_time=LinearTime()
        )
        @test config_linear.flag_time isa LinearTime
    end

    @testset "compute_characteristic_scales" begin
        # Test characteristic scale computation
        μ = 398600.4418
        state = [10000.0, 5000.0, 2000.0, 1.0, 2.0, 3.0]

        DU, TU = AstroCoords.compute_characteristic_scales(state, μ)
        r_mag = sqrt(10000.0^2 + 5000.0^2 + 2000.0^2)
        @test DU ≈ r_mag atol=1e-10
        @test TU ≈ sqrt(r_mag^3 / μ) atol=1e-6

        # Test with zero position (edge case)
        state_zero = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
        DU_zero, TU_zero = AstroCoords.compute_characteristic_scales(state_zero, μ)
        @test DU_zero == 0.0
        @test TU_zero == 0.0
    end

    @testset "compute_initial_phi" begin
        # Test initial phi computation for elliptic orbit (E < 0)
        μ = 398600.4418
        a = 8000.0
        e = 0.1
        # Create a simple circular-ish orbit state
        state = [a, 0.0, 0.0, 0.0, sqrt(μ/a), 0.0]
        config = RegularizedCoordinateConfig(state, μ)

        ϕ₀ = AstroCoords.compute_initial_phi(state, μ, config)
        @test ϕ₀ isa Number
        @test isfinite(ϕ₀)

        # Test with hyperbolic orbit (E > 0)
        v_esc = sqrt(2*μ/a)
        state_hyp = [a, 0.0, 0.0, 0.0, 1.5*v_esc, 0.0]
        config_hyp = RegularizedCoordinateConfig(state_hyp, μ)
        ϕ₀_hyp = AstroCoords.compute_initial_phi(state_hyp, μ, config_hyp)
        @test ϕ₀_hyp isa Number
        @test isfinite(ϕ₀_hyp)

        # Test with perturbing potential
        config_pert = RegularizedCoordinateConfig(state, μ, W=100.0)
        ϕ₀_pert = AstroCoords.compute_initial_phi(state, μ, config_pert)
        @test ϕ₀_pert isa Number
        @test isfinite(ϕ₀_pert)
    end

    @testset "EDromo with RegularizedCoordinateConfig" begin
        # Test EDromo transformations using config
        μ = 398600.4418
        state = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        config = RegularizedCoordinateConfig(state, μ)
        ϕ = AstroCoords.compute_initial_phi(state, μ, config)

        # Test Cartesian to EDromo
        cart = Cartesian(state...)
        edromo = CartesianToEDromo(cart, μ, ϕ, config)
        @test edromo isa EDromo
        @test length(params(edromo)) == 8

        # Test EDromo to Cartesian (round-trip)
        cart_back = EDromoToCartesian(edromo, μ, ϕ, config)
        @test cart_back.x ≈ cart.x atol=1e-6
        @test cart_back.y ≈ cart.y atol=1e-6
        @test cart_back.z ≈ cart.z atol=1e-6
        @test cart_back.ẋ ≈ cart.ẋ atol=1e-6
        @test cart_back.ẏ ≈ cart.ẏ atol=1e-6
        @test cart_back.ż ≈ cart.ż atol=1e-6

        # Test with different flag_time types
        config_linear = RegularizedCoordinateConfig(state, μ, flag_time=LinearTime())
        edromo_linear = CartesianToEDromo(cart, μ, ϕ, config_linear)
        @test edromo_linear isa EDromo

        config_constant = RegularizedCoordinateConfig(state, μ, flag_time=ConstantTime())
        edromo_constant = CartesianToEDromo(cart, μ, ϕ, config_constant)
        @test edromo_constant isa EDromo
    end

    @testset "Kustaanheimo-Stiefel with RegularizedCoordinateConfig" begin
        # Test KS transformations using config
        μ = 398600.4418
        state = [7000.0, 1000.0, 500.0, 0.5, 7.0, 0.2]
        config = RegularizedCoordinateConfig(state, μ, flag_time=LinearTime())

        # Test Cartesian to KS
        cart = Cartesian(state...)
        ks = CartesianToKustaanheimoStiefel(cart, μ, config)
        @test ks isa KustaanheimoStiefel
        @test length(params(ks)) == 10

        # Test KS to Cartesian (round-trip)
        cart_back = KustaanheimoStiefelToCartesian(ks, μ, config)
        @test cart_back.x ≈ cart.x atol=1e-6 rtol=1e-8
        @test cart_back.y ≈ cart.y atol=1e-6 rtol=1e-8
        @test cart_back.z ≈ cart.z atol=1e-6 rtol=1e-8
        @test cart_back.ẋ ≈ cart.ẋ atol=1e-6 rtol=1e-8
        @test cart_back.ẏ ≈ cart.ẏ atol=1e-6 rtol=1e-8
        @test cart_back.ż ≈ cart.ż atol=1e-6 rtol=1e-8
    end

    @testset "Stiefel-Scheifele with RegularizedCoordinateConfig" begin
        # Test SS transformations using config
        μ = 398600.4418
        state = [6500.0, -1000.0, 800.0, 1.0, 6.8, -0.5]
        config = RegularizedCoordinateConfig(state, μ, flag_time=LinearTime())
        ϕ = AstroCoords.compute_initial_phi(state, μ, config)

        # Test Cartesian to SS
        cart = Cartesian(state...)
        ss = CartesianToStiefelScheifele(cart, μ, ϕ, config)
        @test ss isa StiefelScheifele
        @test length(params(ss)) == 10

        # Test SS to Cartesian (round-trip)
        cart_back = StiefelScheifeleToCartesian(ss, μ, ϕ, config)
        @test cart_back.x ≈ cart.x atol=1e-6 rtol=1e-8
        @test cart_back.y ≈ cart.y atol=1e-6 rtol=1e-8
        @test cart_back.z ≈ cart.z atol=1e-6 rtol=1e-8
        @test cart_back.ẋ ≈ cart.ẋ atol=1e-6 rtol=1e-8
        @test cart_back.ẏ ≈ cart.ẏ atol=1e-6 rtol=1e-8
        @test cart_back.ż ≈ cart.ż atol=1e-6 rtol=1e-8
    end
end
