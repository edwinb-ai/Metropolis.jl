mutable struct NVT{V<:Real} <: Ensemble
    δr::V
    accept::V
end

NVT() = NVT(0.5, 0.4)
NVT(x::V) where {V<:Real} = NVT(V(0.5), x)

function _mcmove!(syst::System, uij, opts::EnsembleOptions{T,E}) where {E<:NVT,T}
    (box_size, cutoff) = _box_information(syst.box)
    # Compute the current energy
    uold = _squared_energy(syst.xpos, uij, box_size, cutoff, syst.npart)
    (posold, rng_part) = _choose_move!(
        syst.xpos, syst.rng, opts.ensemble.δr, syst.npartrange
    )
    syst.xpos[rng_part] -= @. box_size * round(syst.xpos[rng_part] / box_size)
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

function _mcmove!(
    syst::System, uij, opts::EnsembleOptions{T,E}, ccache::CellCache; parallel=false
) where {E<:NVT,T}
    # Compute the current energy
    uold = map_pairwise!(uij, 0.0, syst.box, ccache.cell; parallel=parallel)
    (posold, rng_part) = _choose_move!(
        syst.xpos, syst.rng, opts.ensemble.δr, syst.npartrange
    )
    # Update cell lists
    ccache.cell = UpdateCellList!(
        syst.xpos, syst.box, ccache.cell, ccache.aux; parallel=parallel
    )
    # Compute the energy now
    unew = map_pairwise!(uij, 0.0, syst.box, ccache.cell; parallel=parallel)
    Δener = unew - uold

    mcbool = _mcnvt!(Δener, syst.temperature, syst.rng)
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

function _mcnvt!(Δener, temperature, rng)
    if Δener < 0.0 || (rand(rng) < exp(-Δener / temperature))
        # accept += oneunit(accept)
        return true
    else
        return false
    end
end

function _choose_move!(positions, rng, δr, npartrange)
    rng_part = rand(rng, npartrange)
    posold = copy(positions[rng_part])
    # Move that particle
    postrial = 0.5 .- rand!(positions[rng_part])
    positions[rng_part] = @. posold + δr * postrial

    return posold, rng_part
end

function _squared_energy(xpos, uij, box_size, cutoff, npart)
    energy = 0.0

    @inbounds for i in 1:(npart - 1)
        for j in (i + 1):npart
            Δij = xpos[j] - xpos[i]

            # Periodic boundaries
            @. Δij -= box_size * round(Δij / box_size)

            # Compute distance
            Δpos = norm(Δij)^2

            if Δpos < cutoff
                energy = uij(0.0, 0.0, 0, 0, Δpos, energy)
            end
        end
    end

    return energy
end
