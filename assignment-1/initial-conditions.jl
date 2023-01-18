
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

function still_grid_rand(x1, x2, y1, y2, nx, ny, r, m)
    xr = range(x1, x2; length=nx)
    yr = range(y1, y2; length=ny)

    n = nx * ny

    ps = zeros(n, 2)
    vs = zeros(n, 2)
    rs = fill(r, n)
    ms = fill(m, n)

    for (i, x) in enumerate(xr), (j, y) in enumerate(yr)
        ps[begin+ny*(i-1)+(j-1), 1] = x + randn() * r * 0.01
        ps[begin+ny*(i-1)+(j-1), 2] = y + randn() * r * 0.01
    end

    ps, vs, rs, ms
end

# k = packing density.
# k ∈ [0.0, 0.9]
function still_hex_grid(x0, y0, nx, ny, k, d, m)
    r = d * √(√3 / 2π * k)

    d_between = d - 2r

    ps = zeros(2nx * ny, 2)
    vs = zeros(2nx * ny, 2)
    rs = fill(r, 2nx * ny)
    ms = fill(m, 2nx * ny)

    y = y0

    for j in 1:ny
        x = x0
        for i in 1:nx
            ps[begin+2nx*(j-1)+2(i-1), 1] = x + randn() * d_between * 0.1
            ps[begin+2nx*(j-1)+2(i-1), 2] = y + randn() * d_between * 0.1

            ps[begin+2nx*(j-1)+2(i-1)+1, 1] = x + d / 2 +
                                              randn() * d_between * 0.1
            ps[begin+2nx*(j-1)+2(i-1)+1, 2] = y + d * √3 / 2 +
                                              randn() * d_between * 0.1
            x += d
        end

        y += d * √3
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

    ps, vs, rs, ms = still_grid_rand(x1, x2, y1, y2, nx, ny, r, m)

    px = 0.5
    py = p_r + 0.01

    θ = randn() * 0.1 + π / 2
    pvx = p_v * cos(θ)
    pvy = p_v * sin(θ)

    p_area = π * p_r^2
    p_m = p_area * p_density

    [ps; px py], [vs; pvx pvy], [rs; p_r], [ms; p_m]
end

function crater_setup_hex(nx, k, density, p_r, p_density, p_v)
    k′ = √(√3 / 2π * k)

    d = 1 / ((nx - 1 / 2) + 2 - 2k′)
    r = k′ * d

    x = d - r

    ny = round(Int, (1 / 2d - 1) / √3 + 1 / 2)
    y = 1 - ((ny - 1 / 2) * √3 + 1) * d + r

    area = π * r^2
    m = area * density

    ps, vs, rs, ms = still_hex_grid(x, y, nx, ny, k, d, m)

    px = 0.5
    py = p_r + 0.01

    θ = randn() * 0.1 + π / 2
    pvx = p_v * cos(θ)
    pvy = p_v * sin(θ)

    p_area = π * p_r^2
    p_m = p_area * p_density

    [ps; px py], [vs; pvx pvy], [rs; p_r], [ms; p_m]
end

function simple()
    [
        0.5 0.25
        0.5 0.75
    ],
    [
        0.0 1.0
        0.0 0.0
    ],
    [0.03, 0.03],
    [1.0, 1.0]
end
