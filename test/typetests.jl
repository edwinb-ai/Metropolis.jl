@testset "System" begin
    ρ = 0.5
    kT = 1.0
    parts = 10
    cutoff = 1.2

    # Without specifying the cutoff
    s = System(ρ, kT, parts, cutoff)
    @test length(s.xpos) == parts
    @test length(s.xpos[1]) == 3

    # Check that the position vectors are different among them
    @test s.xpos[1] != s.xpos[2] && s.xpos[fld(parts, 2)] != s.xpos[parts]
end

@testset "Simulations" begin
    ρ = 0.5
    kT = 1.0
    parts = 10
    cutoff = 1.2
    s = System(ρ, kT, parts, cutoff)
    sim = Simulation(s, NVT(0.35), LennardJones())

    @test sim isa Simulation
end
