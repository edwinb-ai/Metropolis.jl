"""
    Vec3D{T} <: FieldVector{3,T}

A custom type for mutable static arrays. It completely resembles a 3D point
object in some geometric space.

# Fields
- `x::T`: first field of `Vec3D`
- `y::T`: second field of `Vec3D`
- `z::T`: third field of `Vec3D`
"""
mutable struct Vec3D{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

"Define `similar_type` to preserve array type through array operations."
StaticArrays.similar_type(::Type{<:Vec3D}, ::Type{T}, s::Size{(3,)}) where {T} = Vec3D{T}

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
mutable struct System{B,T,VT,N,R}
    xpos::VT
    density::T
    temperature::T
    box::B
    rng::Random.AbstractRNG
    npart::N
    npartrange::R
end

function System(
    density::T, temp::T, particles::I, cutoff::T; dims=3, lcell=2
) where {T<:Real,I<:Int}
    box_size = cbrt(T(particles) / density)
    box = CellListMap.Box(fill(box_size, dims), cutoff; lcell=lcell)
    rng = Xorshifts.Xoroshiro128Plus()
    xpos = _initialize_positions(box_size, rng, particles)
    syst = System(xpos, density, temp, box, rng, particles, collect(1:particles))

    return syst
end

function _initialize_positions(box_size::T, rng, particles) where {T<:Real}
    range = (zero(T), box_size)
    xpos = [random_vec(Vec3D{T}, range; rng=rng) for _ in 1:particles]

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
