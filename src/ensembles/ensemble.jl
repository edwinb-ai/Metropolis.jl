"""
    Ensemble

A supertype for all kinds of ensembles.
"""
abstract type Ensemble end

mutable struct EnsembleOptions{T,E}
    ensemble::E
    nattempt::T
    naccept::T
end

EnsembleOptions(e::Ensemble) = EnsembleOptions(e, 0, 1)

function mcmove! end
