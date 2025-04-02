using DifferentiationInterface
using Mooncake
using StaticArraysCore

state = 10.0 * randn(6)
state2 = SVector{6}(state)

function test(x::AbstractArray)
    y = SVector{6}(2 * x - x .^ 2)

    return y
end

value_and_jacobian(test, AutoMooncake(; config=nothing), state)[2]
value_and_jacobian(test, AutoMooncake(; config=nothing), state2)

function test2(x::AbstractArray)
    y = Array(2 * x - x .^ 2)

    return y
end

value_and_jacobian(test2, AutoMooncake(; config=nothing), state)
value_and_jacobian(test2, AutoMooncake(; config=nothing), state2)
