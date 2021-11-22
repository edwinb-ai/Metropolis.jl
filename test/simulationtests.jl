function setupsim()
    ρ = 0.78
    kT = 0.85
    particles = 500
    cutoff = 2.5
    s = System(ρ, kT, particles, cutoff; lcell=1)
    return Simulation(s, NVT(0.25, 0.4), LennardJones())
end

@testset "Parallel Simulations" begin
    sim = setupsim()
    @test_nowarn simulate!(sim; steps=10_000, parallel=true)
end

@testset "Serial Simulations" begin
    sim = setupsim()
    @test_nowarn simulate!(sim; steps=10_000, parallel=false)
end
