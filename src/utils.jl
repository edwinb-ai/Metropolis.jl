function random_vec(::Type{VT}, range; rng=Random.GLOBAL_RNG) where {VT}
    dim = length(VT)
    T = eltype(VT)
    (lb, up) = range
    p = VT(lb .+ rand(rng, T, dim) .* (up - lb))

    return p
end

function adjust!(opts::EnsembleOptions)
    if opts.nattempt % opts.naccept == 0
        ratio = opts.naccept / opts.nattempt
        if ratio > opts.ensemble.accept
            opts.ensemble.δr *= 1.05
        else
            opts.ensemble.δr *= 0.95
        end
    end

    return nothing
end
