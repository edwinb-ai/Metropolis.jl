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

# Simulation loop
function _simulation_loop!(system::System, uij, fij, opts, steps; kwargs...)
    for istep in 1:steps
        opts.nattempt += 1
        uener = mcmove!(system, uij, fij, opts, kwargs[:cl]; parallel=kwargs[:parallel])

        if istep % kwargs[:ishow] == 0
            @show uener / system.npart
            @show opts.naccept / opts.nattempt
        end

        adjust!(opts)
    end

    return nothing
end

function _simulate_cells!(sim::Simulation; steps=10_000, parallel=false, ishow=10_000)
    @unpack system, ensemble, potential = sim

    # Obtain the energy function for the interaction potential
    (uij, fij, cl) = _setup_simulation(sim; parallel=parallel)

    # Create the ensemble options
    opts = EnsembleOptions(ensemble)

    _simulation_loop!(system, uij, fij, opts, steps; cl=cl, parallel=parallel, ishow=ishow)

    return nothing
end
