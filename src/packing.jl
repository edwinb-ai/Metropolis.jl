function u_pack(x, box::Box, cl::CellList)
    cl = UpdateCellList!(x, box, cl; parallel=false)
    u = map_pairwise!(
        (x, y, i, j, d2, u) -> begin
            u += (sqrt(d2) - box.cutoff)^2 # objective function
            return u
        end,
        0.0,
        box,
        cl;
        parallel=false,
    )
    return u
end

function gradient_descent!(x::Vector{T}, f, g!; tol=1e-3, maxtrial=500) where {T}
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
