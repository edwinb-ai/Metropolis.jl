"""
    Potential

Supertype for all kinds of interaction potentials.
"""
abstract type Potential end

"""
    Continuous <: Potential

Type for all potentials that implement an potential energy
and a force function. For example, the Lennard-Jones potential.
"""
abstract type Continuous <: Potential end

"""
    Discrete <: Potential

Type for all potentials that only implement a potential energy, due
to its discontinuous nature. For example, the hard-sphere potential.
"""
abstract type Discrete <: Potential end

# Generic definitions for the interface
energy(args...) = nothing
force(args...) = nothing
