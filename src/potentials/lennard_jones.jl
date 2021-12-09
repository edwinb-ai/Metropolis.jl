struct LennardJones <: Continuous
    energy::Function
    force::Function

    function LennardJones()
        pot(x, y, i, j, d2, u) = u += 4.0 * ((1.0 / d2^6) - (1.0 / d2^3))

        function forces(x, y, i, j, d2, f)
            r = x - y
            dudr = -24.0 * r * ((1.0 / d2^7) - (1.0 / d2^4))
            @inbounds f[i] = f[i] + dudr
            @inbounds f[j] = f[j] - dudr

            return f
        end

        return new(pot, forces)
    end
end
