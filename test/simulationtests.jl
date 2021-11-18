@testset "Simulations" begin
    ρ = 0.75
    kT = 1.0
    parts = 1000
    cutoff = 2.5
    s = System(ρ, kT, parts, cutoff)
    sim = Simulation(s, NVT(0.01, 0.35), LennardJones())
    @code_warntype simulate!(sim; steps=100_000)
end
