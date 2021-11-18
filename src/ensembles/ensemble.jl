"""
    Ensemble

A supertype for all kinds of ensembles.
"""
abstract type Ensemble end

mutable struct EnsembleOptions{T,E}
    ensemble::E
    nattempt::T
    naccept::T
    # EnsembleOptions{E,V}(e::E, at::V, acc::V) where {E<:Ensemble, V<:Int} = new(e, at, acc)
end

EnsembleOptions(e::Ensemble) = EnsembleOptions(e, 0, 0)
