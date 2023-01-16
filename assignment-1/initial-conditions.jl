
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

# r = k * d
function crater_setup(nx, k, density)
    r = 1 / (2nx + (nx + 1) / k)
    d = r / k

    x1 = d + r
    x2 = 1 - x1

    ny = nx ÷ 2

    h = 2ny * r + (ny + 1) * d

    y1 = (1 - h) + x1
    y2 = x2

    still_grid(x1, x2, y1, y2, nx, ny, r, 1.0)
end
