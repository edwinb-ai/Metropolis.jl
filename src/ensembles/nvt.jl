mutable struct NVT{V<:Real} <: Ensemble
    δr::V
    accept::V
end

NVT(x::V) where {V<:Real} = NVT(V(0.5), V(x))

function mcmove!(
    syst::System, uij, fij, opts::EnsembleOptions{T,E}; parallel=false
) where {E<:NVT,T<:Real}
    (box_size, cutoff) = _box_information(syst.box)
    # Compute the current energy
    uold = _compute_energy(syst.xpos, uij, box_size, cutoff, syst.npart)
    (posold, rng_part) = _choose_move!(syst.xpos, syst.npart, syst.rng, opts.ensemble.δr)
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
    syst::System, uij, fij, opts::EnsembleOptions{T,E}, cl; parallel=false
) where {E<:NVT,T<:Real}
    # Compute the current energy
    uold = map_pairwise!(uij, zero(T), syst.box, cl; parallel=parallel)
    (posold, rng_part) = _choose_move!(syst.xpos, syst.npart, syst.rng, opts.ensemble.δr)
    # Update cell lists
    cl = UpdateCellList!(syst.xpos, syst.box, cl; parallel=parallel)
    # Compute the energy now
    unew = map_pairwise!(uij, zero(T), syst.box, cl; parallel=parallel)
    Δener = unew - uold

    mcbool = _mcnvt!(unew, uold, Δener, syst.temperature, opts.naccept; rng=syst.rng)
    if mcbool
        uold += Δener
        opts.naccept += 1
    else
        syst.xpos[rng_part] = posold
        # Update cell lists
        cl = UpdateCellList!(syst.xpos, syst.box, cl; parallel=parallel)
    end

    return uold
end

function _mcnvt!(unew, uold, dener, temperature, accept; rng=Random.GLOBAL_RNG)
    if unew < uold
        if rand(rng) < exp(-dener / temperature)
            return true
        end
    else
        return false
    end
end

function _choose_move!(positions, npart, rng, δr)
    rng_part = rand(rng, 1:(npart))
    posold = positions[rng_part]
    vec_size = length(posold)
    # Move that particle
    half_disp = Vec3D(0.5 .- rand(rng, vec_size))
    positions[rng_part] = posold .+ δr .* half_disp

    return posold, rng_part
end

function _compute_energy(xpos, uij, box_size, cutoff, npart)
    energy = 0.0

    @inbounds for i in 1:(npart - 1)
        for j in (i + 1):npart
            Δij = xpos[i] - xpos[j]

            # Periodic boundaries
            pbc_ij = @. Δij - box_size * round(Δij / box_size)

            # Compute distance
            Δpos = LinearAlgebra.norm(pbc_ij)

            if Δpos < cutoff
                energy = uij(0.0, 0.0, 0, 0, Δpos, energy)
            end
        end
    end

    return energy
end
