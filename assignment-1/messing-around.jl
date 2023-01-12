using SIMD

function find_intersection_time(p1x, p1y, p2x, p2y, v1x, v1y, v2x, v2y, r1, r2)
    Δpx = p2x - p1x
    Δpy = p2y - p1y

    Δvx = v2x - v1x
    Δvy = v2y - v1y

    a = Δvx^2 + Δvy^2
    b = Δpx * Δvx + Δpy * Δvy
    c = Δpx^2 + Δpy^2 - (r1 + r2)^2

    @fastmath (-b - sqrt(b^2 - a * c)) / a
end

function find_first_collision_single(ps, vs, rs, i)
    min_ind = 0
    min_t = NaN

    @inbounds for j in (i+1):size(ps, 1)
        t = @inline find_intersection_time(
            ps[i, 1], ps[i, 2], ps[j, 1], ps[j, 2],
            vs[i, 1], vs[i, 2], vs[j, 1], vs[j, 2],
            rs[i], rs[j]
        )

        if t > 0 && !(min_t < t)
            min_ind = j
            min_t = t
        end
    end

    min_ind, min_t
end

function range_chunks(r, n)
    rv = @inbounds r[1:n:end-(n-1)]
    rr = @inbounds r[length(rv)*n+1:end]
    rv, rr
end

function rising_vec(V::Type{Vec{N,T}}) where {N,T<:Integer}
    V(((0:N-1)...,))
end

function find_first_collision_single_simd(ps, vs, rs, i,
    V::Type{Vec{N,Float64}}) where {N}
    VI = Vec{N,Int}

    j_range = (i+1):size(ps, 1)

    jv, jr = @inline range_chunks(j_range, N)

    min_ind = VI(0)
    min_t = V(NaN)

    @inbounds for j in jv
        jv = VecRange{N}(j)
        t = @inline find_intersection_time(
            V(ps[i, 1]), V(ps[i, 2]), ps[jv, 1], ps[jv, 2],
            V(vs[i, 1]), V(vs[i, 2]), vs[jv, 1], vs[jv, 2],
            V(rs[i]), rs[jv]
        )

        mask = (t > 0) & !(min_t < t)
        min_ind = vifelse(mask, VI(j) + rising_vec(VI), min_ind)
        min_t = vifelse(mask, t, min_t)
    end

    min_ind_single = 0
    min_t_single = NaN

    @inbounds for lane in 1:N
        j = min_ind[lane]
        t = min_t[lane]

        if t > 0 && !(min_t_single < t)
            min_ind_single = j
            min_t_single = t
        end
    end

    @inbounds for j in jr
        t = @inline find_intersection_time(
            ps[i, 1], ps[i, 2], ps[j, 1], ps[j, 2],
            vs[i, 1], vs[i, 2], vs[j, 1], vs[j, 2],
            rs[i], rs[j]
        )

        if t > 0 && !(min_t_single < t)
            min_ind_single = j
            min_t_single = t
        end
    end

    min_ind_single, min_t_single
end

function find_first_collision(ps, vs, rs)
    min_inds = (0, 0)
    min_t = NaN

    @inbounds for i in axes(ps, 1)
        j, t = @inline find_first_collision_single_simd(
            ps, vs, rs, i,
            Vec{8,Float64}
        )
        # j, t = @inline find_first_collision_single(ps, vs, rs, i)

        if t > 0 && !(min_t < t)
            min_inds = (i, j)
            min_t = t
        end
    end

    min_inds, min_t
end

function find_first_collision_threaded(ps, vs, rs)
    min_inds = [(0, 0) for _ in 1:Threads.nthreads()]
    min_t = [NaN for _ in 1:Threads.nthreads()]

    Threads.@threads for i in axes(ps, 1)
        j, t = @inline find_first_collision_single(ps, vs, rs, i)

        thid = Threads.threadid()

        if t > 0 && !(min_t[thid] < t)
            min_inds[thid] = (i, j)
            min_t[thid] = t
        end
    end

    min_inds_single = (0, 0)
    min_t_single = NaN

    for (inds, t) in zip(min_inds, min_t)
        if t > 0 && !(min_t_single < t)
            min_inds_single = inds
            min_t_single = t
        end
    end

    min_inds_single, min_t_single
end

function setup_still_square_lattice(n, r, spacing)
    ps = zeros(n^2, 2)
    vs = zeros(n^2, 2)
    rs = fill(r, n^2)

    for i in 0:n-1, j in 0:n-1
        ps[begin+n*i+j, 1] = Float64(i) * spacing
        ps[begin+n*i+j, 2] = Float64(j) * spacing
    end

    ps, vs, rs
end

function setup_square_bowling(n, r, spacing, px, py, vx, vy)
    ps, vs, rs = setup_still_square_lattice(n, r, spacing)

    [ps; px py], [vs; vx vy], [rs; r]
end
