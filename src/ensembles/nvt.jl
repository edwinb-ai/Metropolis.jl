mutable struct NVT{V<:Real} <: Ensemble
    δr::V
    accept::V
end

NVT(x::V) where {V<:Real} = NVT(V(0.5), V(x))

function mcmove!(syst::System, uij, fij, opts::EnsembleOptions{T,E}) where {E<:NVT,T<:Real}
    (box_size, cutoff) = _box_information(syst.box)
    # Compute the current energy
    uold = _compute_energy(syst.xpos, uij, box_size, cutoff, syst.npart)
    (posold, rng_part) = _choose_move!(
        syst.xpos, syst.rng, opts.ensemble.δr, syst.npartrange
    )
    # Compute the energy now
    unew = _compute_energy(syst.xpos, uij, box_size, cutoff, syst.npart)
    Δener = unew - uold

    mcbool = _mcnvt!(unew, uold, Δener, syst.temperature, opts.naccept; rng=syst.rng)
    if mcbool
        uold += Δener
        opts.naccept += 1
    else
        syst.xpos[rng_part] = posold
    end

    return uold
end

function mcmove!(
    syst::System, uij, fij, opts::EnsembleOptions{T,E}, cl::CellList; parallel=false
) where {E<:NVT,T<:Int}
    # Compute the current energy
    uold = map_pairwise!(uij, 0.0, syst.box, cl; parallel=parallel)
    (posold, rng_part) = _choose_move!(
        syst.xpos, syst.rng, opts.ensemble.δr, syst.npartrange
    )
    # Update cell lists
    cl = UpdateCellList!(syst.xpos, syst.box, cl; parallel=parallel)
    # Compute the energy now
    unew = map_pairwise!(uij, 0.0, syst.box, cl; parallel=parallel)
    Δener = unew - uold

    mcbool = _mcnvt!(unew, uold, Δener, syst.temperature, opts.naccept; rng=syst.rng)
    if mcbool
        uold += Δener
        opts.naccept += one(T)
    else
        syst.xpos[rng_part] = posold
        # Update cell lists
        cl = UpdateCellList!(syst.xpos, syst.box, cl; parallel=parallel)
    end

    return uold
end

function _mcnvt!(unew, uold, dener, temperature, accept; rng=Random.GLOBAL_RNG)
    if unew < uold || (rand(rng) < exp(-dener / temperature))
        accept += oneunit(accept)
        return true
    else
        return false
    end
end

function _choose_move!(positions, rng, δr, npartrange)
    rng_part = rand(rng, npartrange)
    posold = copy(positions[rng_part])
    # Move that particle
    half_disp = Vec3D(0.5 .- rand(rng, 3))
    positions[rng_part] = @. posold + δr * half_disp

    return posold, rng_part
end

function _compute_energy(xpos, uij, box_size, cutoff, npart)
    energy = 0.0

    @inbounds for i in 1:(npart - 1)
        for j in (i + 1):npart
            Δij = xpos[i] - xpos[j]

            # Periodic boundaries
            Δij = @. Δij - box_size * fld(Δij, box_size)

            # Compute distance
            Δpos = LinearAlgebra.norm(Δij)^2

            if Δpos < cutoff
                energy = uij(0.0, 0.0, 0, 0, Δpos, energy)
            end
        end
    end

    return energy
end
