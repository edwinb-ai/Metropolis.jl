function random_vec(::Type{VT}, range; rng=Random.GLOBAL_RNG) where {VT}
    dim = length(VT)
    T = eltype(VT)
    (lb, up) = range
    p = VT(lb .+ rand(rng, T, dim) .* (up - lb))

    return p
end
