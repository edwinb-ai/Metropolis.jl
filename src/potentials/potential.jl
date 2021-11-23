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
struct Continuous <: Potential end

"""
    Discrete <: Potential

Type for all potentials that only implement a potential energy, due
to its discontinuous nature. For example, the hard-sphere potential.
"""
struct Discrete <: Potential end
