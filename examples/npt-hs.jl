using Metropolis
using Printf
using JLD2

function run(syst, dispm, simul; pos=nothing)
    # Create the positions vector
    positions = [zeros(SVector{3}) for n in 1:syst.N]
    if isnothing(pos)
        # Initialize the positions as grid
        init!(positions, syst)
    else
        positions = copy(pos)
    end

    save_file = @sprintf "npt_%.3f_%.3f.csv" syst.ρ syst.P
    positions_file = @sprintf "positions_%.3f_%.2f.jld2" syst.ρ syst.P

    # * Main loop
    # Equilibration steps
    move(positions, syst, simul.eq, dispm; volume=true)
    JLD2.@save joinpath("data", positions_file) positions


    # Sampling steps
    move(
        positions,
        syst,
        simul.sample,
        dispm;
        volume=true,
        filename=joinpath("results", save_file),
    )
end

(syst, dispm, simul, posfile) = parse_toml("npt.toml")

run(syst, dispm, simul; pos=posfile)
