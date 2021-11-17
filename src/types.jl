"""
    Ensemble

A supertype for all kinds of ensembles.
"""
abstract type Ensemble end

"""
    System{V<:Real}

Holds all relevant information from the simulation.
"""
mutable struct System{V, VT}
    density::V
    temperature::V
    xpos::VT
end
