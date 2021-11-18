mutable struct NVT{V<:Real} <: Ensemble
    δr::V
    accept::V
end

NVT(x::V) where {V<:Real} = NVT(V(0.5), V(x))

function mcmove!(
    syst::System, uij, fij, opts::EnsembleOptions{T,E}, cl; parallel=false
) where {E<:NVT,T<:Real}
    @unpack ensemble, nattempt, naccept = opts

    # Compute the current energy
    uold = map_pairwise!(uij, 0.0, syst.box, cl; parallel=parallel)
    # Choose a random particle
    rng_part = rand(syst.rng, 1:(syst.npart))
    # Save this particle's position
    posold = syst.xpos[rng_part]
    vec_size = length(posold)
    # Move that particle
    half_disp = SVector{vec_size,eltype(posold)}(0.5 .- rand(syst.rng, vec_size))
    new_pos = @. posold + ensemble.δr * half_disp
    syst.xpos[rng_part] = new_pos
    # Update cell lists
    cl = UpdateCellList!(syst.xpos, syst.box, cl; parallel=parallel)
    # Compute the energy now
    unew = map_pairwise!(uij, 0.0, syst.box, cl; parallel=parallel)
    Δener = unew - uold

    if unew < uold
        if rand(syst.rng) < exp(-Δener / syst.temperature)
            uold += Δener
            naccept += 1
        end
    else
        syst.xpos[rng_part] = posold
        # Update cell lists
        cl = UpdateCellList!(syst.xpos, syst.box, cl; parallel=parallel)
    end

    @pack! opts = naccept

    return uold, cl
end
