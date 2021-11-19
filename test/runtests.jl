using Metropolis
using Test

@testset "Metropolis.jl" begin
    include("test_interface.jl")
    include("typetests.jl")
    include("simulationtests.jl")
end
