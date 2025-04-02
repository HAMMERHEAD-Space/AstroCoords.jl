for set in _COORDINATE_SETS
    @test length(check_allocs(set, (Vector{Number},))) == 0
end

for set1 in _COORDINATE_SETS
    for set2 in _COORDINATE_SETS
        @test length(check_allocs(set1, (set2{Float64}, Float64))) == 0
    end
end

for set in _COORDINATE_SETS
    @test length(check_allocs(meanMotion, (set{Float64}, Float64))) == 0
    @test length(check_allocs(orbitalPeriod, (set{Float64}, Float64))) == 0
    @test length(check_allocs(orbitalNRG, (set{Float64}, Float64))) == 0
    @test length(check_allocs(angularMomentumVector, (set{Float64}, Float64))) == 0
    @test length(check_allocs(angularMomentumQuantity, (set{Float64}, Float64))) == 0
end
