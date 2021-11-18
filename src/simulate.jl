function simulate! end

function simulate!(sim::Simulation; steps=10_000, parallel=false, ishow=10_000)
    @unpack system, ensemble, potential = sim

    # Obtain the energy function for the interaction potential
    (uij, fij) = interactions(potential)

    # Build initial cell lists
    cl = CellList(copy(system.xpos), system.box; parallel=parallel)

    # Create the ensemble options
    opts = EnsembleOptions(ensemble)

    # Create a good enough configuration
    packsystem!(system, cl, uij)

    # * Simulation loop
    for istep in 1:steps
        opts.nattempt += 1
        (upot, cl) = mcmove!(system, uij, fij, opts, cl)

        if istep % ishow == 0
            @show upot / system.npart
            @show opts.naccept / opts.nattempt
        end

        adjust!(opts)
    end

    return nothing
end
