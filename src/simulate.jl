function simulate! end

function simulate!(sim::Simulation; pairwise=:cells, kwargs...)
    if pairwise == :cells
        _simulate_cells!(sim; kwargs...)
    else
        _simulate_squared!(sim; kwargs...)
    end

    return nothing
end

function _setup_simulation(sim::Simulation; kwargs...)
    # Obtain the energy function for the interaction potential
    (uij, fij) = interactions(sim.potential)

    # Build initial cell lists
    cl = CellListMap.CellList(sim.system.xpos, sim.system.box; parallel=kwargs[:parallel])

    # Create a good enough configuration
    packsystem!(sim.system, cl, uij)

    return uij, fij, cl
end

function _simulate_cells!(sim::Simulation; steps=10_000, parallel=false, ishow=10_000)
    # Obtain the energy function for the interaction potential
    (uij, _, cl) = _setup_simulation(sim; parallel=parallel)

    # Create the ensemble options
    opts = EnsembleOptions(sim.ensemble)

    # for istep in 1:steps
    @showprogress for istep in 1:steps
        opts.nattempt += oneunit(opts.nattempt)
        (uener, cl) = mcmove!(sim.system, uij, opts, cl; parallel=parallel)

        if istep % ishow == 0
            @show uener / sim.system.npart
            @show opts.naccept / opts.nattempt
        end

        adjust!(opts)
    end

    return nothing
end

function _simulate_squared!(sim::Simulation; steps=10_000, parallel=false, ishow=10_000)
    # Obtain the energy function for the interaction potential
    (uij, fij, _) = _setup_simulation(sim; parallel=parallel)

    # Create the ensemble options
    opts = EnsembleOptions(sim.ensemble)

    @showprogress for istep in 1:steps
        opts.nattempt += 1
        uener = mcmove!(sim.system, uij, fij, opts)

        if istep % ishow == 0
            @show uener / sim.system.npart
            @show opts.naccept / opts.nattempt
        end

        adjust!(opts)
    end

    return nothing
end
