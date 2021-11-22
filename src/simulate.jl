function simulate! end

function simulate!(sim::Simulation; pairwise=:cells, kwargs...)
    if pairwise == :cells
        _simulate_cells!(sim; kwargs...)
    else
        _simulate_squared!(sim; kwargs...)
    end

    return nothing
end

function _setup_simulation(system::System, potential::Potential; kwargs...)
    # Obtain the energy function for the interaction potential
    (uij, fij) = interactions(potential)

    # Build initial cell lists
    cl = CellListMap.CellList(system.xpos, system.box; parallel=kwargs[:parallel])

    # Create a good enough configuration
    packsystem!(system, cl, uij)

    return uij, fij, cl
end

function _simulate_cells!(sim::Simulation; steps=10_000, parallel=false, ishow=10_000)
    @unpack system, ensemble, potential = sim
    # Obtain the energy function for the interaction potential
    (uij, _, cl) = _setup_simulation(system, potential; parallel=parallel)

    # Create the ensemble options
    opts = EnsembleOptions(ensemble)
    aux = CellListMap.AuxThreaded(cl)

    for istep in 1:steps
    # @showprogress for istep in 1:steps
        opts.nattempt += 1
        (uener, cl) = mcmove!(system, uij, opts, cl, aux; parallel=parallel)

        if istep % ishow == 0
            @show uener / system.npart
            @show opts.naccept / opts.nattempt
        end

        adjust!(opts)
    end

    return nothing
end

function _simulate_squared!(sim::Simulation; steps=10_000, parallel=false, ishow=10_000)
    @unpack system, ensemble, potential = sim
    # Obtain the energy function for the interaction potential
    (uij, _, _) = _setup_simulation(system, potential; parallel=parallel)

    # Create the ensemble options
    opts = EnsembleOptions(sim.ensemble)

    for istep in 1:steps
        opts.nattempt += 1
        uener = mcmove!(sim.system, uij, opts)

        if istep % ishow == 0
            @show uener / sim.system.npart
            @show opts.naccept / opts.nattempt
        end

        adjust!(opts)
    end

    return nothing
end
