function simulate! end

function simulate!(sim::Simulation; steps=10_000, parallel=false)
    @unpack system, ensemble, potential = sim
    # @unpack xpos, density, temperature, box, rng, npart = system

    # Obtain the energy function for the interaction potential
    uij = potential_energy(potential)
    fij = forces(potential)

    # Build initial cell lists
    cl = CellList(copy(system.xpos), system.box.box; parallel=parallel)

    # Create the ensemble options
    opts = EnsembleOptions(ensemble)

    uenergy = map_pairwise!(uij, 0.0, system.box.box, cl) / system.npart
    println("initial energy $(uenergy)")

    # * Simulation loop
    for istep in 1:steps
        opts.nattempt += 1
        (upot, cl) = mcmove!(system, uij, fij, opts, cl)
        uenergy += upot

        if istep % 1_000 == 0
            @show uenergy / system.npart
            @show upot / system.npart
            @show opts.naccept / opts.nattempt
        end

        adjust!(opts)
    end

    return nothing
end
