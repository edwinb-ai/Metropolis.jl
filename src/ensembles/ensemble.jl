"""
    Ensemble

A supertype for all kinds of ensembles.
"""
abstract type Ensemble end

struct EnsembleOptions{T,V,E}
    ensemble::E
    nattempt::T
    naccept::T
end

EnsembleOptions(e::Ensemble) = EnsembleOptions(e, 0, 0)
