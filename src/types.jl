"""
    Ensemble

A supertype for all kinds of ensembles.
"""
abstract type Ensemble end

"""
    System{V<:Real}

Holds all relevant information from the simulation.
"""
mutable struct System{UnitCellType,N,T,M,VT}
    xpos::VT
    density::T
    temperature::T
    box::CellListMap.Box{UnitCellType,N,T,M}
    cutoff::T
    rng::Random.AbstractRNG
end
