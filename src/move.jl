function is_overlap(pos::AbstractArray, syst::System)
    rc = syst.σ^2
    @inbounds for i = 1:syst.N - 1
        for j = (i + 1):syst.N
            Δij = pos[i] - pos[j]

            # Periodic boundaries
            pbc_ij = @. Δij - syst.L * round(Δij / syst.L)

            # Compute distance
            Δpos2 = norm(pbc_ij)^2

            # If there is overlap, return immediately
            if Δpos2 < rc
                return true
            end
        end
    end

    # If no overlaps encountered, return false for the full configuration
    return false
end

function count_overlap(pos::AbstractArray, syst::System)
    rc = syst.σ^2
    @inbounds for i = 1:syst.N - 1
        for j = (i + 1):syst.N
            Δij = pos[i] - pos[j]

            # Periodic boundaries
            pbc_ij = @. Δij - syst.L * round(Δij / syst.L)

            # Compute distance
            Δpos2 = norm(pbc_ij)^2

            if Δpos2 < rc
                syst.n_overlap += 1
            end
        end
    end

    return nothing
end

function move(
    positions::AbstractArray,
    syst::System,
    cycles::Int,
    disp::Displacements;
    volume::Bool=false,
    filename=nothing,
)
    # Run values
    attempts = 0
    volatt = 0
    accepted = 0
    volaccpt = 0
    total_ρ = 0
    ratio = 0
    volratio = 0
    samples = 0

    # Equilibration steps
    @showprogress for i = 1:cycles
        # * Attempt to move the particles around
        attempts += 1
        mc = metropolis!(positions, syst, disp)
        accepted += mc
        ratio = accepted / attempts

        # * Attempt to change the volume in the system
        if volume
            if i % Int(2 * syst.N) == 0
                volatt += 1
                mc = mcvolume!(positions, syst, disp)
                volaccpt += mc
                volratio = volaccpt / volatt
            end
        end

        # * Save to file
        if i % 1_00 == 0
            samples += 1

            # * Update the total density of the system
            if volume
                total_ρ += syst.N / syst.vol
            end

            if !isnothing(filename)
                open(filename, "a") do io
                    if volume
                        filerho = total_ρ / samples
                    else
                        filerho = syst.ρ
                    end
                    println(io, "$(i),$(filerho)")
                end
            end
        end
    end
    # Save previous values
    if volume
        syst.ρ = total_ρ / samples
    end

    # Adjust results as averages
    if volume
        total_ρ /= samples
    end

    # * Show results
    println("Ratio of acceptance: $(ratio)")
    if volume
        println("Density: $(total_ρ)")
        println("Volume ratio: $(volratio)")
    end
end
