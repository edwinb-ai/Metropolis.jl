"""
    Ensemble

A supertype for all kinds of ensembles.
"""
abstract type Ensemble end

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
    box::CellListMap.Box{UnitCellType,N,T,M}
    cutoff::T
    rng::Random.AbstractRNG
    npart::Int
end

function System(density::T, temp::T, particles::Int; dims=3) where {T<:Real}
    box_size = particles / density
    cutoff = box_size / 2.0
    box = CellListMap.Box(fill(box_size, dims), cutoff; lcell=2)
    rng = Xorshifts.Xoroshiro128Plus()
    xpos = fill(random_vec(SVector{3,Float64}, (zero(T), box_size); rng=rng), particles)
    syst = System(xpos, density, temp, box, cutoff, rng, particles)

    return syst
end

function System(density::T, temp::T, particles::Int, cutoff::T; dims=3) where {T<:Real}
    box_size = particles / density
    box = CellListMap.Box(fill(box_size, dims), cutoff; lcell=2)
    rng = Xorshifts.Xoroshiro128Plus()
    xpos = fill(random_vec(SVector{3,Float64}, (zero(T), box_size); rng=rng), particles)
    syst = System(xpos, density, temp, box, cutoff, rng, particles)

    return syst
end
