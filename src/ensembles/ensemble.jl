"""
    Ensemble

A supertype for all kinds of ensembles.
"""
abstract type Ensemble end

"""
    EnsembleOptions{T,E}

A type to hold necessary information for each kind of `Ensemble`.

# Fields
- `ensemble::E`: the ensemble that will specify this type; normally a sub-type of
    `Ensemble`.
- `nattempt::T`: the number of attempts that each ensemble has been called.
- `naccept::T`: the number of accepted attempts for each ensemble called.
"""
mutable struct EnsembleOptions{T,E}
    ensemble::E
    nattempt::T
    naccept::T
end

EnsembleOptions(e::Ensemble) = EnsembleOptions(e, 0, 1)
