"""
    Potential

Supertype for all kinds of interaction potentials.
"""
abstract type Potential end

"""
    System{V<:Real}

Holds all relevant information for the simulation system.

# Fields
- `xpos::VT`: stores the array for the positions of the particles.
- `density::T`: the density of the system.
- `temperature::T`: the temperature of the system.
- `box::CellListMap.Box{UnitCellType,N,T,M}`: simulation box that will get updated as the
    simulation runs.
- `cutoff::T`: cutoff radius for the molecular interactions
- `rng::Random.AbstractRNG`: holds the RNG object for the system
- `npart::Int`: total number of particles in the system
"""
mutable struct System{UnitCellType,N,T,M,VT}
    xpos::VT
    density::T
    temperature::T
    sbox::CellListMap.Box{UnitCellType,N,T,M}
    rng::Random.AbstractRNG
    npart::Int
end

function System(
    density::T, temp::T, particles::Int, cutoff::T; dims=3, random_init=true
) where {T<:Real}
    box_size = cbrt(particles / density)
    box = CellListMap.Box(fill(box_size, dims), cutoff; lcell=2)
    rng = Xorshifts.Xoroshiro128Plus()
    xpos = initialize_positions(box, rng, particles; dims=dims, random_init=random_init)
    syst = System(xpos, density, temp, box, rng, particles)

    return syst
end

function initialize_positions(box_size, rng, particles; dims=3, random_init=true)
    if random_init
        range = (zero(T), box_size)
        xpos = [random_vec(SVector{dims,Float64}, range; rng=rng) for _ in 1:particles]
    else
        xpos = [zeros(SVector{dims,Float64}) for _ in 1:particles]
        square_lattice!(xpos, particles, box_size)
    end

    return xpos
end

"""
    Simulation

A type to hold all relevant information for the simulation.

# Fields
- `system::S`: ideally, a `System` type that has all the information about the positions
    and physical information on the system.
- `ensemble::E`: this is the _statistical ensemble_ used for doing the simulation, ideally
    of type `Ensemble`.
- `potential::P`: the interaction potential between the particles, ideally an object of
"""
struct Simulation{S,E,P}
    system::S
    ensemble::E
    potential::P
end
