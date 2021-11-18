function random_vec(::Type{VT}, range; rng=Random.GLOBAL_RNG) where {VT}
    dim = length(VT)
    T = eltype(VT)
    (lb, up) = range
    p = VT(lb .+ rand(rng, T, dim) .* (up - lb))

    return p
end

function square_lattice!(positions, boxl, npart)
    n3 = 2
    ix = 0
    iy = 0
    iz = 0

    # Find the lowest perfect cube, n3, greater than or equal to the
    # number of particles
    while n3^3 < npart
        n3 += 1
    end

    for i in axes(positions, 1)
        new_pos = Vec3D(
            (ix + 0.5) * boxl / n3, (iy + 0.5) * boxl / n3, (iz + 0.5) * boxl / n3
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

    return nothing
end

function adjust!(opts::EnsembleOptions)
    @unpack ensemble, nattempt, naccept = opts

    if nattempt % naccept == 0
        ratio = naccept / nattempt
        if ratio > ensemble.accept
            ensemble.δr *= 1.05
        else
            ensemble.δr *= 0.95
        end
    end

    @pack! opts = ensemble, nattempt, naccept

    return nothing
end
