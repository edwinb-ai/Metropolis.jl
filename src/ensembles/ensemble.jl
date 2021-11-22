"""
    Ensemble

A supertype for all kinds of ensembles.
"""
abstract type Ensemble end

mutable struct EnsembleOptions{T,E,V}
    ensemble::E
    nattempt::T
    naccept::T
    movecache::V
end

EnsembleOptions(e::Ensemble) = EnsembleOptions(e, 0, 1, zeros(MVector{3, Float64}))

function mcmove! end
