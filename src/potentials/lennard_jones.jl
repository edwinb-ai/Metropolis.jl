struct LennardJones <: Continuous
    energy::Function
    force::Function

    function LennardJones()
        pot(d2, σ) = 4.0 * ((σ^12 / d2^6) - (σ^6 / d2^3))
        forces(d2, σ, r) = -24.0 * r * ((σ^12 / d2^7) - (σ^6 / d2^4))

        return new(pot, forces)
    end
end

# Implement the generic interfaces
function potential_energy(lj::LennardJones)
    lj_energy(d2, σ, u) = u += lj.energy(d2, σ)

    return lj_energy
end

function forces(lj::LennardJones)
    function lj_forces(x, y, i, j, d2, σ, f)
        r = x - y
        dudr = lj.forces(d2, σ, r)
        @inbounds f[i] = f[i] + dudr
        @inbounds f[j] = f[j] - dudr

        return f
    end

    return lj_forces
end
