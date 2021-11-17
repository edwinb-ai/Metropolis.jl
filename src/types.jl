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
    npart::Int
end

function System{UnitCellType,N,T,M,VT}(
    density::T, temp::T, particles::Int; dims=3
) where {UnitCellType,N,T,M,VT}
    box_size = particles / density
    # volume = cbrt(box_size)
    cutoff = box_size / 2.0
    box = CellListMap.Box(fill(box_size, dims), cutoff; lcell=2)
    rng = Xorshifts.Xoroshiro128Plus()
    xpos = fill(random_vec(SVector{3,Float64}, (zero(T), box_size)), particles)
    syst = System(xpos, density, temp, box, cutoff, rng, particles)

    return syst
end
