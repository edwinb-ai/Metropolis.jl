using Metropolis
using Test

@testset "Metropolis.jl" begin
    lj = LennardJones()
    @test potential_energy(lj) isa Function
    @test forces(lj) isa Function
end
