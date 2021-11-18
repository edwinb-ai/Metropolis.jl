function simulate! end

function simulate!(sim::Simulation; steps=10_000, parallel=false)
    @unpack system, ensemble, potential = sim

    # Obtain the energy function for the interaction potential
    uij = potential_energy(potential)
    fij = forces(potential)

    # Build initial cell lists
    x = copy(system.xpos)
    cl = CellList(x, system.box; parallel=parallel)

    # Create the ensemble options
    opts = EnsembleOptions(ensemble)

    uenergy = map_pairwise!(uij, 0.0, system.box, cl)

    # * Simulation loop
    for istep in 1:steps
        (upot, cl) = mcmove!(system, ensemble, uij, fij, opts, cl)
        uenergy += upot

        if istep % 1_000 == 0
            @show uenergy / system.npart
        end
    end

    return nothing
end
