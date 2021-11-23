function simulate! end

function simulate!(sim::Simulation; pairwise=:cells, kwargs...)
    if pairwise == :cells
        _simulate_cells!(sim; kwargs...)
    else
        _simulate_squared!(sim; kwargs...)
    end

    return nothing
end

function _setup_simulation!(system::System, potential::Potential; kwargs...)
    # Build initial cell lists
    cl = CellListMap.CellList(system.xpos, system.box; parallel=kwargs[:parallel])
    cell_cache = _build_cache(cl)

    # Create a good enough configuration
    packsystem!(system, cell_cache.cell, potential.energy)

    return cell_cache
end

function _simulate_cells!(sim::Simulation; steps=10_000, parallel=false, ishow=10_000)
    @unpack system, ensemble, potential = sim
    # Obtain the energy function for the interaction potential
    cell_cache = _setup_simulation!(system, potential; parallel=parallel)

    # Create the ensemble options
    opts = EnsembleOptions(ensemble)

    @showprogress for istep in 1:steps
        opts.nattempt += 1
        uener = mcmove!(system, potential.energy, opts, cell_cache; parallel=parallel)

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
    _ = _setup_simulation!(system, potential; parallel=parallel)

    # Create the ensemble options
    opts = EnsembleOptions(ensemble)

    for istep in 1:steps
        opts.nattempt += 1
        uener = mcmove!(system, potential.energy, opts)

        if istep % ishow == 0
            @show uener / system.npart
            @show opts.naccept / opts.nattempt
        end

        adjust!(opts)
    end

    return nothing
end
