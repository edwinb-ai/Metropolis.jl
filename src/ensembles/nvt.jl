mutable struct NVT{V<:Real} <: Ensemble
    δr::V
    accept::V
end

NVT(x::V) where {V<:Real} = NVT(V(0.5), V(x))

function mcmove!(
    syst::System, ens::E, uij, fij, opts::EnsembleOptions{T,E}, cl
) where {E<:NVT,T<:Real}
    # Compute the current energy
    uold = map_pairwise!(uij, 0.0, syst.box, cl)
    # Choose a random particle
    rng_part = rand(syst.rng, 1:(syst.npart))
    # Save this particle's position
    posold = syst.xpos[rng_part]
    vec_size = length(posold)
    # Move that particle
    half_disp = SVector{vec_size,eltype(posold)}(0.5 .- rand(syst.rng, vec_size))
    new_pos = @. posold + ens.δr * half_disp
    syst.xpos[rng_part] = new_pos
    # Update cell lists
    cl = UpdateCellList!(syst.xpos, syst.box, cl)
    # Compute the energy now
    unew = map_pairwise!(uij, 0.0, syst.box, cl)
    Δener = unew - uold

    opts.nattempt += 1
    if unew < uold
        if rand(syst.rng) < exp(-Δener / syst.temperature)
            uold = Δener
            opts.naccept += 1
        end
    else
        syst.xpos[rng_part] = posold
        # Update cell lists
        cl = UpdateCellList!(syst.xpos, syst.box, cl)
    end

    return uold, cl
end
