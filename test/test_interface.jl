@testset "Lennard Jones" begin
    lj = LennardJones()
    @test potential_energy(lj) isa Function
    @test forces(lj) isa Function
end

@testset "Hard Sphere" begin
    struct HardSphere <: Discrete
        energy::Function

        HardSphere() = new((d2, σ) -> Inf)
    end

    function potential_energy(hs::HardSphere)
        hs_energy(d2, σ, u) = u += hs.energy(d2, σ)

        return hs_energy
    end

    hs = HardSphere()
    @test potential_energy(hs) isa Function
end
