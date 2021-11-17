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
function potential_energy end
function forces end
