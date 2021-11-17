@testset "System" begin
    ρ = 0.5
    kT = 1.0
    parts = 1000

    # Without specifying the cutoff
    s = System(ρ, kT, parts)
    @test length(s.xpos) == parts
    @test length(s.xpos[1]) == 3

    # Specifying the cutoff
    cutoff = 3.0
    s = System(ρ, kT, parts, cutoff)
    @test s.cutoff == cutoff
end
