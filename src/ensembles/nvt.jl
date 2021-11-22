mutable struct NVT{V<:Real} <: Ensemble
    δr::V
    accept::V
end

NVT(x::V) where {V<:Real} = NVT(V(0.5), V(x))

function mcmove!(syst::System, uij, opts::EnsembleOptions{T,E}) where {E<:NVT,T}
    (box_size, cutoff) = _box_information(syst.box)
    # Compute the current energy
    uold = _squared_energy(syst.xpos, uij, box_size, cutoff, syst.npart)
    (posold, rng_part) = _choose_move!(
        syst.xpos, syst.rng, opts.ensemble.δr, syst.npartrange, opts.movecache
    )
    # Compute the energy now
    unew = _squared_energy(syst.xpos, uij, box_size, cutoff, syst.npart)
    Δener = unew - uold

    mcbool = _mcnvt!(unew, uold, Δener, syst.temperature, opts.naccept, syst.rng)
    if mcbool
        uold += Δener
        opts.naccept += oneunit(T)
    else
        syst.xpos[rng_part] = posold
    end

    return uold
end

function mcmove!(
    syst::System, uij, opts::EnsembleOptions{T,E}, ccache::CellCache; parallel=false
) where {E<:NVT,T}
    # Compute the current energy
    uold = map_pairwise!(uij, 0.0, syst.box, ccache.cell; parallel=parallel)
    (posold, rng_part) = _choose_move!(
        syst.xpos, syst.rng, opts.ensemble.δr, syst.npartrange, opts.movecache
    )
    # Update cell lists
    ccache.cell = UpdateCellList!(
        syst.xpos, syst.box, ccache.cell, ccache.aux; parallel=parallel
    )
    # Compute the energy now
    unew = map_pairwise!(uij, 0.0, syst.box, ccache.cell; parallel=parallel)
    Δener = unew - uold

    mcbool = _mcnvt!(unew, uold, Δener, syst.temperature, opts.naccept, syst.rng)
    if mcbool
        uold += Δener
        opts.naccept += oneunit(T)
    else
        syst.xpos[rng_part] = posold
        # Update cell lists
        ccache.cell = UpdateCellList!(
            syst.xpos, syst.box, ccache.cell, ccache.aux; parallel=parallel
        )
    end

    return uold
end

function _mcnvt!(unew, uold, dener, temperature, accept, rng)
    if unew < uold || (rand(rng) < exp(-dener / temperature))
        accept += oneunit(accept)
        return true
    else
        return false
    end
end

function _choose_move!(positions, rng, δr, npartrange, movevec)
    rng_part = rand(rng, npartrange)
    posold = copy(positions[rng_part])
    # Move that particle
    half_disp = 0.5 .- rand!(rng, movevec)
    positions[rng_part] = @. posold + δr * half_disp

    return posold, rng_part
end

function _squared_energy(xpos, uij, box_size, cutoff, npart)
    energy = 0.0
    rc = cutoff^2

    @inbounds for i in 1:(npart - 1)
        for j in (i + 1):npart
            Δij = xpos[i] - xpos[j]

            # Periodic boundaries
            Δij = @. Δij - box_size * round(Δij / box_size)

            # Compute distance
            Δpos = norm(Δij)^2

            if Δpos < rc
                energy = uij(0.0, 0.0, 0, 0, Δpos, energy)
            end
        end
    end

    return energy
end
