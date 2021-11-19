@testset "Simulations" begin
    ρ = 0.7727
    kT = 0.85
    particles = 500
    cutoff = 2.5
    s = System(ρ, kT, particles, cutoff; lcell=1)
    sim = Simulation(s, NVT(0.25, 0.4), LennardJones())
    @test_nowarn simulate!(sim; steps=10_000)
end
