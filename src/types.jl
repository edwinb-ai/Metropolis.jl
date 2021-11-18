"""
    Potential

Supertype for all kinds of interaction potentials.
"""
abstract type Potential end

struct Vec3D{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

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
mutable struct System{B,T,VT,I}
    xpos::VT
    density::T
    temperature::T
    box::B
    rng::Random.AbstractRNG
    npart::I
end

function System(
    density::T, temp::T, particles::I, cutoff::T; dims=3, random_init=true, lcell=2
) where {T<:Real,I<:Int}
    box_size = cbrt(T(particles) / density)
    # box = create_box(box_size, dims, cutoff; lcell=lcell)
    box = CellListMap.Box(fill(box_size, dims), cutoff; lcell=lcell)
    rng = Xorshifts.Xoroshiro128Plus()
    xpos = initialize_positions(box_size, rng, particles; random_init=random_init)
    syst = System(xpos, density, temp, box, rng, particles)

    return syst
end

function create_box(box_size::T, dims::I, cutoff::T; lcell=2) where {T,I}
    box = CellListMap.Box(fill(box_size, dims), cutoff; lcell=lcell)

    return box
end

function initialize_positions(box_size::T, rng, particles; random_init=true) where {T<:Real}
    xpos = [zeros(Vec3D{T}) for _ in 1:particles]
    if random_init
        range = (zero(typeof(box_size)), box_size)
        xpos = [random_vec(Vec3D{T}, range; rng=rng) for _ in 1:particles]
    else
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
