using DataStructures

function find_intersection_time(p1x, p1y, p2x, p2y, v1x, v1y, v2x, v2y, r1, r2)
    Δpx = p2x - p1x
    Δpy = p2y - p1y

    Δvx = v2x - v1x
    Δvy = v2y - v1y

    @fastmath begin
        a = Δvx^2 + Δvy^2
        b = Δpx * Δvx + Δpy * Δvy
        d = r1 + r2
        c = Δpx^2 + (Δpy + d) * (Δpy - d)

        (-b - sqrt(b^2 - a * c)) / a
    end
end

@inline function sorttuple((i, j),)
    i <= j ? (i, j) : (j, i)
end

function delete_occurences!(pq, i, n)
    # Delete all occurences of the index i in the pq
    # Time should be O(n) + O(n_o * log(n))
    # where n_o is number of occurences of the index i
    for j in -3:n
        ij = sorttuple((i, j))
        if haskey(pq, ij)
            delete!(pq, ij)
        end
    end
end

# indices for walls:
#  0: right wall    (+x)
# -1: left wall     (-x)
# -2: bottom wall   (+y)
# -3: top wall      (-y)
function add_wall_collisions!(pq, i, t, ps, vs, rs)
    x = ps[i, 1]
    y = ps[i, 2]

    vx = vs[i, 1]
    vy = vs[i, 2]

    r = rs[i]

    Δtx = ((vx > 0 ? 1 - r : r) - x) / vx
    Δty = ((vy > 0 ? 1 - r : r) - y) / vy

    j = 0
    Δt = Inf

    if isfinite(Δtx)
        j = (vx > 0) - 1
        Δt = Δtx
    end

    if isfinite(Δty) && Δty < Δt
        j = (vy > 0) - 3
        Δt = Δty
    end

    ij = (j, i)
    if isfinite(Δt)
        pq[ij] = t + Δt
    elseif haskey(pq, ij)
        delete!(pq, ij)
    end
end

function add_disk_collisions!(pq, i, t, ps, vs, rs)
    n = size(ps, 1)

    for j in 1:n
        Δt = find_intersection_time(
            ps[i, 1], ps[i, 2], ps[j, 1], ps[j, 2],
            vs[i, 1], vs[i, 2], vs[j, 1], vs[j, 2],
            rs[i], rs[j]
        )

        ij = sorttuple((i, j))
        if isfinite(Δt) && Δt > 0
            pq[ij] = t + Δt
        elseif haskey(pq, ij)
            delete!(pq, ij)
        end
    end
end

function update_collisions!(pq, i, t, ps, vs, rs)
    if i >= 1
        # delete_occurences!(pq, i, size(ps, 1))
        add_wall_collisions!(pq, i, t, ps, vs, rs)
        add_disk_collisions!(pq, i, t, ps, vs, rs)
    end
end

function init_collisions(ps, vs, rs)
    pq = PriorityQueue{NTuple{2,Int},Float64}()

    for i in axes(ps, 1)
        update_collisions!(pq, i, 0.0, ps, vs, rs)
    end

    pq
end
