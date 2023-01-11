using SIMD

function find_intersection_time(p1, p2, v1, v2, r1, r2)
    Δp = p2 .- p1
    Δv = v2 .- v1

    a = Δv[1]^2 + Δv[2]^2
    b = 2 * (Δp[1] * Δv[1] + Δp[2] * Δv[2])
    c = Δp[1]^2 + Δp[2]^2 - (r1 + r2)^2

    @fastmath 0.5 * (-b - sqrt(b^2 - 4a * c)) / a
end
