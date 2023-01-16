
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
function crater_setup(nx, k, density, p_r, p_density, p_v)
    r = 1 / (2nx + (nx + 1) / k)
    d = r / k

    x1 = d + r
    x2 = 1 - x1

    ny = nx ÷ 2

    h = 2ny * r + (ny + 1) * d

    y1 = (1 - h) + x1
    y2 = x2

    area = π * r^2
    m = area * density

    ps, vs, rs, ms = still_grid(x1, x2, y1, y2, nx, ny, r, m)

    px = 0.5
    py = p_r + 0.01

    θ = randn() * 0.1 + π / 2
    pvx = p_v * cos(θ)
    pvy = p_v * sin(θ)

    p_area = π * p_r^2
    p_m = p_area * p_density

    [ps; px py], [vs; pvx pvy], [rs; p_r], [ms; p_m]
end
