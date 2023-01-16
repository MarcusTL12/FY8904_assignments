
function still_grid(x1, x2, y1, y2, nx, ny, r, m)
    xr = range(x1, x2; length=nx)
    yr = range(y1, y2; length=ny)

    n = nx * ny

    ps = zeros(n, 2)
    vs = zeros(n, 2)
    rs = fill(r, n)
    ms = fill(m, n)

    for (i, x) in enumerate(xr), (j, y) in enumerate(yr)
        ps[begin+ny*(i-1)+(j-1), 1] = x
        ps[begin+ny*(i-1)+(j-1), 2] = y
    end

    ps, vs, rs, ms
end
