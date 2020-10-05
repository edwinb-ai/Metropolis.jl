using Random: AbstractRNG

mutable struct Displacements{V <: Real}
    δV::V
    δr::V
end

mutable struct Simulation{V <: Int}
    eq::V
    sample::V
end

mutable struct System{V <: Real}
    N::Int
    L::V
    σ::V
    ρ::V
    vol::V
    press::V
    β::V
    P::V
    n_overlap::Int
    rng::AbstractRNG
end
