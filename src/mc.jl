function metropolis!(pos::AbstractArray, syst::System, disp::Displacements)
    # Choose a random particle
    rng_part = rand(syst.rng, 1:syst.N)
    # Save this particle's position
    posold = pos[rng_part]
    # Move that particle
    half_disp = 0.5 .- rand(syst.rng, 3)
    new_pos = @. pos[rng_part] + disp.δr * half_disp
    pos[rng_part] = new_pos
    # Apply periodic boundary conditions
    @inbounds for (i, p) in enumerate(pos)
        new_pos = @. p - syst.L * round(p / syst.L)
        pos[i] = new_pos
    end

    # Check for overlap
    # overlap = overlap_tree(pos, syst)
    overlap = is_overlap(pos, syst)
    if overlap # Reject the movement of particles
        pos[rng_part] = posold
        return 0
    else # Accept it
        return 1
    end
end

function mcvolume!(pos::AbstractArray, syst::System, disp::Displacements)
    # Initial volume
    vol_old = syst.vol
    # Obtain the natural log
    lnVold = log(vol_old)
    # Compute the new volume
    lnvolnew = lnVold + disp.δV * (rand(syst.rng) - 0.5)
    vol_new = exp(lnvolnew)
    # Compute the new box length
    Lnew = ∛vol_new

    # Adjust the particles to new box
    adjust = Lnew / syst.L
    syst.L = Lnew
    for (i, p) in enumerate(pos)
        new_pos = p * adjust
        pos[i] = new_pos
    end

    # Compute new density
    ρold = syst.ρ
    syst.ρ = syst.N / vol_new

    # Check for overlaps
    # overlap = overlap_tree(pos, syst)
    overlap = is_overlap(pos, syst)
    if overlap # Reject the movement of particles
        for (i, p) in enumerate(pos)
            new_pos = p / adjust
            pos[i] = new_pos
        end
        syst.ρ = ρold
        syst.L = ∛vol_old
        syst.vol = vol_old
        return 0
    # We use the Metropolis criteria
    else
        # Compute the energy difference
        ΔH = syst.P * (vol_new - vol_old) - (syst.N + 1) * log(vol_new / vol_old)
        # Apply Metropolis' criteria
        if rand(syst.rng) <= exp(-ΔH)
            # println("Volume accepted!")
            syst.vol = vol_new
            return 1
        else
            # println("Volume rejected!")
            # Re-scale the positions
            for (i, p) in enumerate(pos)
                new_pos = p / adjust
                pos[i] = new_pos
            end
            syst.ρ = ρold
            syst.L = ∛vol_old
            syst.vol = vol_old
            return 0
        end
    end
end
