@testset "Types" begin
    ρ = 0.5
    kT = 1.0
    parts = 1000
    s = System(ρ, kT, parts)

    @test length(s.xpos) == parts
    @test length(s.xpos[1]) == 3
end
