struct LennardJones <: Continuous
    energy::Function
    force::Function

    function LennardJones()
        pot(d2) = 4.0 * ((1.0 / d2^6) - (1.0 / d2^3))
        forces(d2, r) = -24.0 * r * ((1.0 / d2^7) - (1.0 / d2^4))

        return new(pot, forces)
    end
end

# Implement the generic interfaces
function potential_energy(lj::LennardJones)
    lj_energy(x, y, i, j, d2, u) = u = lj.energy(d2)

    return lj_energy
end

function forces(lj::LennardJones)
    function lj_forces(x, y, i, j, d2, f)
        r = x - y
        dudr = lj.forces(d2, r)
        @inbounds f[i] = f[i] + dudr
        @inbounds f[j] = f[j] - dudr

        return f
    end

    return lj_forces
end
