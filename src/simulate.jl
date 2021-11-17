function simulate!(sim::Simulation; steps=10_000, parallel=false)
    @unpack syst, ens, pot = sim

    # Obtain the energy function for the interaction potential
    uij = potential_energy(pot)

    # Build initial cell lists
    x = copy(syst.xpos)
    cl = CellList(x, syst.box; parallel=parallel)

    # Create the ensemble options

    # * Simulation loop
    for istep in 1:steps
        mcmove!(syst, ens, pot)

        if istep % 1_000 == 0
            total_energy = map_pairwise!(uij, 0.0, syst.box, cl; parallel=parallel)
            @show total_energy
        end
    end

    return nothing
end
