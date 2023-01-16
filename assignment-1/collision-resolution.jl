using SIMD

function do_wall_collision!(i, j, vs, ξ)
    if -1 <= i
        vs[j, 1] *= -ξ
    else
        vs[j, 2] *= -ξ
    end
end

# Updates velocities of particle i and j (or only j if i is wall)
function do_collision!(i, j, ps, vs, ms, ξ)
    if i < 1
        do_wall_collision!(i, j, vs, ξ)
        return
    end

    p_i = Vec((ps[i, 1], ps[i, 2]))
    p_j = Vec((ps[j, 1], ps[j, 2]))

    Δp = p_j - p_i

    Rij2 = Δp[1]^2 + Δp[2]^2

    vi = Vec((vs[i, 1], vs[i, 2]))
    vj = Vec((vs[j, 1], vs[j, 2]))

    Δv = vj - vi

    mi = ms[i]
    mj = ms[j]

    common = (1 + ξ) * sum(Δv * Δp) * Δp / ((mi + mj) * Rij2)

    vi_new = vi + mj * common
    vj_new = vj - mi * common

    vs[i, 1] = vi_new[1]
    vs[i, 2] = vi_new[2]

    vs[j, 1] = vj_new[1]
    vs[j, 2] = vj_new[2]
end
