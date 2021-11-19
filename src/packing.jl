function u_pack(x, box::Box, cl::CellList; parallel=false)
    cl = UpdateCellList!(x, box, cl; parallel=parallel)
    u = map_pairwise!(
        (x, y, i, j, d2, u) -> begin
            u += (sqrt(d2) - box.cutoff)^2 # objective function
            return u
        end,
        0.0,
        box,
        cl;
        parallel=parallel,
    )
    return u
end

function gradient_descent!(x::Vector{T}, f, g!; tol=1e-5, maxtrial=50_000) where {T}
    gnorm(x) = maximum(norm(v) for v in x)
    itrial = 0
    step = 1.0
    xtrial = similar(x)
    g = fill!(similar(x), zero(T))
    fx = f(x)
    g = g!(g, x)
    while (gnorm(g) > tol) && (itrial < maxtrial)
        @. xtrial = x - step * g
        ftrial = f(xtrial)
        if ftrial >= fx
            step = step / 2
        else
            x .= xtrial
            fx = ftrial
            g = g!(g, x)
            step = step * 2
        end
        itrial += 1
    end
    return x
end

function fpair_cl(x, y, i, j, d2, f, box::Box)
    Δv = y - x
    d = sqrt(d2)
    fₓ = 2.0 * (d - box.cutoff) * (Δv / d)
    f[i] += fₓ
    f[j] -= fₓ
    return f
end

function forces_cl!(f::Vector{T}, x, box::Box, cl::CellList, fpair::F) where {T,F}
    fill!(f, zero(T))
    cl = UpdateCellList!(x, box, cl; parallel=false)
    map_pairwise!(
        (x, y, i, j, d2, f) -> fpair(x, y, i, j, d2, f, box), f, box, cl; parallel=false
    )
    return f
end

function packpositions(positions, box::Box; parallel=false)
    box_pack = Box(box.unit_cell_max, 0.5)
    cl_pack = CellList(positions, box_pack; parallel=parallel)
    xpos = gradient_descent!(
        copy(positions),
        (x) -> u_pack(x, box_pack, cl_pack),
        (g, x) -> -forces_cl!(g, x, box_pack, cl_pack, fpair_cl),
    )

    return xpos, box_pack, cl_pack
end

function packsystem!(system::System, cl::CellList, uij; parallel=false)
    uenergy = map_pairwise!(uij, 0.0, system.box, cl; parallel=parallel) / system.npart
    println("initial energy $(uenergy)")

    (pack_pos, boxpack, clpack) = packpositions(copy(system.xpos), system.box)
    cl = UpdateCellList!(pack_pos, system.box, cl; parallel=parallel)
    upacked = u_pack(pack_pos, boxpack, clpack; parallel=parallel)
    system.xpos = copy(pack_pos)
    @show upacked

    uenergy = map_pairwise!(uij, 0.0, system.box, cl; parallel=parallel) / system.npart
    println("initial energy $(uenergy)")

    return nothing
end
