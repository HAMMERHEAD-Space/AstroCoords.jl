@testset "Cartesian Constructors" begin
    @testset "Inner constructor" begin
        # Test line 26: Cartesian{T}(...)
        cart_f64 = Cartesian{Float64}(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
        @test cart_f64 isa Cartesian{Float64}
        @test cart_f64.x ≈ 1.0 atol=1e-15
        @test cart_f64.y ≈ 2.0 atol=1e-15
        @test cart_f64.z ≈ 3.0 atol=1e-15
        @test cart_f64.ẋ ≈ 4.0 atol=1e-15
        @test cart_f64.ẏ ≈ 5.0 atol=1e-15
        @test cart_f64.ż ≈ 6.0 atol=1e-15

        cart_f32 = Cartesian{Float32}(1.0f0, 2.0f0, 3.0f0, 4.0f0, 5.0f0, 6.0f0)
        @test cart_f32 isa Cartesian{Float32}
        @test eltype(cart_f32) == Float32
    end

    @testset "Constructor from AbstractVector" begin
        # Test lines 33-35: Cartesian(X::AbstractVector)
        vec = [7.0, 8.0, 9.0, 1.0, 2.0, 3.0]
        cart = Cartesian(vec)
        @test cart isa Cartesian
        @test cart.x ≈ 7.0 atol=1e-15
        @test cart.y ≈ 8.0 atol=1e-15
        @test cart.z ≈ 9.0 atol=1e-15
        @test cart.ẋ ≈ 1.0 atol=1e-15
        @test cart.ẏ ≈ 2.0 atol=1e-15
        @test cart.ż ≈ 3.0 atol=1e-15

        # Test with different vector types
        vec_f32 = Float32[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        cart_f32 = Cartesian(vec_f32)
        @test eltype(cart_f32) == Float32
    end

    @testset "Constructor with mixed types (promoting)" begin
        # Test line 37: type promotion
        cart_mixed = Cartesian(1, 2.0, 3.0f0, 4, 5.0, 6.0f0)
        @test cart_mixed isa Cartesian
        # Should promote to Float64
        @test eltype(cart_mixed) == Float64
        @test cart_mixed.x ≈ 1.0 atol=1e-15
        @test cart_mixed.y ≈ 2.0 atol=1e-15
        @test cart_mixed.z ≈ 3.0 atol=1e-15

        # All integers
        cart_int = Cartesian(1, 2, 3, 4, 5, 6)
        @test eltype(cart_int) == Int
    end

    @testset "Constructor from StaticVector" begin
        # Test line 40: (::Type{C})(g::StaticVector)
        svec = SVector{6}(10.0, 11.0, 12.0, 13.0, 14.0, 15.0)
        cart = Cartesian(svec)
        @test cart isa Cartesian
        @test cart.x ≈ 10.0 atol=1e-15
        @test cart.y ≈ 11.0 atol=1e-15
        @test cart.z ≈ 12.0 atol=1e-15
        @test cart.ẋ ≈ 13.0 atol=1e-15
        @test cart.ẏ ≈ 14.0 atol=1e-15
        @test cart.ż ≈ 15.0 atol=1e-15
    end

    @testset "Base.one(Cartesian)" begin
        # Test lines 43-46: Base.one with different T types
        cart_one_f64 = Base.one(Cartesian; T=Float64)
        @test cart_one_f64 isa Cartesian{Float64}
        @test cart_one_f64.x ≈ 0.0 atol=1e-15
        @test cart_one_f64.y ≈ 0.0 atol=1e-15
        @test cart_one_f64.z ≈ 0.0 atol=1e-15
        @test cart_one_f64.ẋ ≈ 0.0 atol=1e-15
        @test cart_one_f64.ẏ ≈ 0.0 atol=1e-15
        @test cart_one_f64.ż ≈ 0.0 atol=1e-15

        cart_one_f32 = Base.one(Cartesian; T=Float32)
        @test cart_one_f32 isa Cartesian{Float32}
        @test eltype(cart_one_f32) == Float32

        cart_one_int = Base.one(Cartesian; T=Int64)
        @test cart_one_int isa Cartesian{Int64}
        @test eltype(cart_one_int) == Int64
    end

    @testset "Base.getindex for all 6 indices" begin
        # Test lines 48-60: indexing
        cart = Cartesian(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)

        @test cart[1] ≈ 1.0 atol=1e-15  # x
        @test cart[2] ≈ 2.0 atol=1e-15  # y
        @test cart[3] ≈ 3.0 atol=1e-15  # z
        @test cart[4] ≈ 4.0 atol=1e-15  # ẋ
        @test cart[5] ≈ 5.0 atol=1e-15  # ẏ
        @test cart[6] ≈ 6.0 atol=1e-15  # ż

        # Test out of bounds
        @test_throws BoundsError cart[0]
        @test_throws BoundsError cart[7]
        @test_throws BoundsError cart[-1]
    end

    @testset "Type stability" begin
        # Ensure all constructors are type-stable
        vec = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        cart = Cartesian(vec)
        @test @inferred Cartesian(vec) isa Cartesian{Float64}
        @test @inferred Cartesian(1.0, 2.0, 3.0, 4.0, 5.0, 6.0) isa Cartesian{Float64}
        @test @inferred Base.one(Cartesian; T=Float64) isa Cartesian{Float64}
    end

    @testset "params function" begin
        # Test that params returns SVector correctly
        cart = Cartesian(7.0, 8.0, 9.0, 1.0, 2.0, 3.0)
        p = params(cart)
        @test p isa SVector{6,Float64}
        @test p[1] ≈ 7.0 atol=1e-15
        @test p[2] ≈ 8.0 atol=1e-15
        @test p[3] ≈ 9.0 atol=1e-15
        @test p[4] ≈ 1.0 atol=1e-15
        @test p[5] ≈ 2.0 atol=1e-15
        @test p[6] ≈ 3.0 atol=1e-15
    end

    @testset "Property getters: r and v" begin
        # Test lines 71-80: getproperty for :r and :v
        cart = Cartesian(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)

        # Test position vector property :r
        r = cart.r
        @test r isa SVector{3,Float64}
        @test r[1] ≈ 1.0 atol=1e-15  # x
        @test r[2] ≈ 2.0 atol=1e-15  # y
        @test r[3] ≈ 3.0 atol=1e-15  # z

        # Test velocity vector property :v
        v = cart.v
        @test v isa SVector{3,Float64}
        @test v[1] ≈ 4.0 atol=1e-15  # ẋ
        @test v[2] ≈ 5.0 atol=1e-15  # ẏ
        @test v[3] ≈ 6.0 atol=1e-15  # ż

        # Test with different numeric types
        cart_f32 = Cartesian{Float32}(1.0f0, 2.0f0, 3.0f0, 4.0f0, 5.0f0, 6.0f0)
        @test cart_f32.r isa SVector{3,Float32}
        @test cart_f32.v isa SVector{3,Float32}

        # Test that regular fields still work
        @test cart.x ≈ 1.0 atol=1e-15
        @test cart.ẋ ≈ 4.0 atol=1e-15

        # Test type stability
        @test @inferred(cart.r) isa SVector{3,Float64}
        @test @inferred(cart.v) isa SVector{3,Float64}
    end

    @testset "propertynames" begin
        # Test line 82-84: propertynames
        cart = Cartesian(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
        pnames = propertynames(cart)

        @test :x in pnames
        @test :y in pnames
        @test :z in pnames
        @test :ẋ in pnames
        @test :ẏ in pnames
        @test :ż in pnames
        @test :r in pnames
        @test :v in pnames

        # Verify it's the expected tuple
        @test pnames == (:x, :y, :z, :ẋ, :ẏ, :ż, :r, :v)
    end
end
