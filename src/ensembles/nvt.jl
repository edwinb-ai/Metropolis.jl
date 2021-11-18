mutable struct NVT{V<:Real} <: Ensemble
    δr::V
    accept::V
end

NVT(x::V) where {V<:Real} = NVT(V(0.5), V(x))

function _choose_move!(positions, npart, rng, δr)
    rng_part = rand(rng, 1:(npart))
    posold = positions[rng_part]
    vec_size = length(posold)
    # Move that particle
    half_disp = Vec3D(0.5 .- rand(rng, vec_size))
    positions[rng_part] = posold .+ δr .* half_disp

    return posold, rng_part
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

    if unew < uold
        if rand(syst.rng, T) < exp(-Δener / syst.temperature)
            uold += Δener
            opts.naccept += 1
        end
    else
        syst.xpos[rng_part] = posold
        # Update cell lists
        cl = UpdateCellList!(syst.xpos, syst.box, cl; parallel=parallel)
    end

    return uold
end
