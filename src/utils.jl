function init!(positions::AbstractArray, syst::System)
    n3 = 2
    ix = 0
    iy = 0
    iz = 0

    # Find the lowest perfect cube, n3, greater than or equal to the
    # number of particles
    while n3^3 < syst.N
        n3 += 1
    end

    for i in axes(positions, 1)
        new_pos = SVector(
            (ix + 0.5) * syst.L / n3, (iy + 0.5) * syst.L / n3, (iz + 0.5) * syst.L / n3
        )
        positions[i] = new_pos
        ix += 1

        if ix == n3
            ix = 0
            iy += 1
            if iy == n3
                iy = 0
                iz += 1
            end
        end
    end
end

function tofile(x::String, filename::String)
    open(filename, "a") do io
        println(io, x)
    end
end

function parse_toml(s::String)
    result = TOML.parsefile(s)

    syst = System(0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, MersenneTwister())

    for (i, j) in result["system"]
        if i == "particles"
            syst.N = j
        elseif i == "pressure"
            syst.P = j
        elseif i == "density"
            syst.ρ = j
        elseif i == "temperature"
            syst.β == 1 / j
        end
    end

    syst.vol = syst.N / syst.ρ
    syst.L = ∛(syst.vol)

    if "seed" in keys(result["system"])
        syst.rng = Xoroshiro128Plus(result["system"]["seed"])
    else
        syst.rng = Xoroshiro128Plus()
    end

    dispm = Displacements(1.0, 0.5)

    for (i, j) in result["displacements"]
        if i == "position"
            dispm.δr = j
        elseif i == "logvolume"
            dispm.δV = j
        end
    end

    simulation = Simulation(100000, 200000)

    for (i, j) in result["simulation"]
        if i == "equilibration"
            simulation.eq = j
        elseif i == "sampling"
            simulation.sample = j
        end
    end

    if "title" in keys(result)
        println(result["title"])
        println("Initial density: $(syst.ρ)")
        println("Pressure: $(syst.P)")
    end

    positions = zeros(3, syst.N)
    if "positions" in keys(result)
        println("Reading positions from file...")
        JLD2.@load result["positions"]["file"] positions
    else
        positions = nothing
    end

    return (syst, dispm, simulation, positions)
end
