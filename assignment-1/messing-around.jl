using SIMD

function find_intersection_time(p1x, p1y, p2x, p2y, v1x, v1y, v2x, v2y, r1, r2)
    Δpx = p2x - p1x
    Δpy = p2y - p1y

    Δvx = v2x - v1x
    Δvy = v2y - v1y

    @fastmath begin
        a = Δvx^2 + Δvy^2
        b = Δpx * Δvx + Δpy * Δvy
        # c = Δpx^2 + Δpy^2 - (r1 + r2)^2
        d = r1 + r2
        c = Δpx^2 + (Δpy + d) * (Δpy - d)

        (-b - sqrt(b^2 - a * c)) / a
    end
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

function find_first_collision_simd(ps, vs, rs)
    N = 8
    V = Vec{N,Float64}
    VI = Vec{N,Int}

    min_inds = (0, 0)
    min_t = NaN

    i_range = 2:size(ps, 1)

    iv, ir = @inline range_chunks(i_range, N)

    @inbounds for i in iv
        iv = VecRange{N}(i)
        pxi = ps[iv, 1]
        pyi = ps[iv, 2]
        vxi = vs[iv, 1]
        vyi = vs[iv, 2]
        ri = rs[iv]

        min_i_v = VI(0)
        min_j_v = VI(0)
        min_t_v = V(NaN)

        for j in i+N:size(ps, 1)
            t = @inline find_intersection_time(
                pxi, pyi, V(ps[j, 1]), V(ps[j, 2]),
                vxi, vyi, V(vs[j, 1]), V(vs[j, 2]),
                ri, V(rs[j])
            )

            mask = (t > 0) & !(min_t_v < t)
            min_i_v = vifelse(mask, VI(i) + rising_vec(VI), min_i_v)
            min_j_v = vifelse(mask, VI(j), min_j_v)
            min_t_v = vifelse(mask, t, min_t_v)
        end

        min_i = 0
        min_j = 0
        min_t_local = NaN

        for lane in 1:N
            i_lane = min_i_v[lane]
            j_lane = min_j_v[lane]
            t_lane = min_t_v[lane]

            if t_lane > 0 && !(min_t_local < t_lane)
                min_i = i_lane
                min_j = j_lane
                min_t_local = t_lane
            end
        end

        for ii in i:i+N-2, j in ii+1:i+N-1
            t = @inline find_intersection_time(
                ps[ii, 1], ps[ii, 2], ps[j, 1], ps[j, 2],
                vs[ii, 1], vs[ii, 2], vs[j, 1], vs[j, 2],
                rs[ii], rs[j]
            )

            if t > 0 && !(min_t_local < t)
                min_i = ii
                min_j = j
                min_t_local = t
            end
        end

        if min_t_local > 0 && !(min_t < min_t_local)
            min_inds = (min_i, min_j)
            min_t = min_t_local
        end
    end

    @inbounds for i in ir
        j, t = @inline find_first_collision_single(ps, vs, rs, i)

        if t > 0 && !(min_t < t)
            min_inds = (i, j)
            min_t = t
        end
    end

    min_inds, min_t
end

function find_first_collision_threaded(ps, vs, rs)
    nth = Threads.nthreads()

    min_inds = [(0, 0) for _ in 1:nth]
    min_t = [NaN for _ in 1:nth]

    Threads.@threads for id in 1:nth
        for i in id:nth:size(ps, 1)
            j, t = @inline find_first_collision_single_simd(
                ps, vs, rs, i,
                Vec{8,Float64}
            )
            # j, t = @inline find_first_collision_single(ps, vs, rs, i)

            if t > 0 && !(min_t[id] < t)
                min_inds[id] = (i, j)
                min_t[id] = t
            end
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
