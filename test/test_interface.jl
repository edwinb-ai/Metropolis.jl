@testset "Lennard Jones" begin
    lj = LennardJones()
    @test lj.energy isa Function
    @test lj.force isa Function
end

@testset "Hard Sphere" begin
    struct HardSphere <: Discrete
        energy::Function

        HardSphere() = new((x, y, i, j, d2, u) -> u += Inf)
    end

    hs = HardSphere()
    @test hs.energy isa Function
end
