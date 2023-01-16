using Images
using SIMD

include("simd-utils.jl")

function show_disks(ps, rs, resolution)
    img = zeros(resolution, resolution)

    xr = range(0.0, 1.0; length=resolution)
    yr = range(0.0, 1.0; length=resolution)

    for (j, x) in enumerate(xr), (i, y) in enumerate(yr)
        for ((px, py), r) in zip(eachrow(ps), rs)
            d = √((x - px)^2 + (y - py)^2)
            Δd = (d - r) * resolution
            if Δd <= -0.5
                img[i, j] = 1.0
                break
            elseif Δd <= 0.5
                img[i, j] = max(0.5 - Δd, img[i, j])
            end
        end
    end

    Gray.(img)
end

function show_disks_simd(ps, rs, resolution)
    N = 8
    V = Vec{N,Float64}

    img = zeros(resolution, resolution)

    xr = range(0.0, 1.0; length=resolution)
    yr = range(0.0, 1.0; length=resolution)

    krv, krr = range_chunks(1:size(ps, 1), N)

    @inbounds for (j, x) in enumerate(xr), (i, y) in enumerate(yr)
        pixel_v = V(img[i, j])

        for k in krv
            kv = VecRange{N}(k)
            pxv = ps[kv, 1]
            pyv = ps[kv, 2]
            rv = rs[kv]

            Δdv = (√((x - pxv)^2 + (y - pyv)^2) - rv) * resolution

            pixel_v = vifelse(Δdv <= -0.5, 1.0, pixel_v)
            pixel_v = vifelse((-0.5 < Δdv) & (Δdv <= 0.5),
                max(0.5 - Δdv, pixel_v), pixel_v)
        end

        img[i, j] = maximum(pixel_v)

        for k in krr
            px = ps[k, 1]
            py = ps[k, 2]
            r = rs[k]

            Δd = (√((x - px)^2 + (y - py)^2) - r) * resolution

            if Δd <= -0.5
                img[i, j] = 1.0
                break
            elseif Δd <= 0.5
                img[i, j] = max(0.5 - Δd, img[i, j])
            end
        end
    end

    Gray.(img)
end

function show_disks_smart(ps, rs, resolution)
    img = zeros(resolution, resolution)

    xr = range(0.0, 1.0, resolution)
    yr = range(0.0, 1.0, resolution)

    @inbounds for ((px, py), r) in zip(eachrow(ps), rs)
        x1 = px - r
        x2 = px + r
        j1 = clamp(floor(Int, x1 * resolution), 1, resolution)
        j2 = clamp(ceil(Int, x2 * resolution), 1, resolution)

        y1 = py - r
        y2 = py + r
        i1 = clamp(floor(Int, y1 * resolution), 1, resolution)
        i2 = clamp(ceil(Int, y2 * resolution), 1, resolution)

        for j in j1:j2, i in i1:i2
            x = xr[j]
            y = yr[i]

            d = √((x - px)^2 + (y - py)^2)
            Δd = (d - r) * resolution
            if Δd <= -0.5
                img[i, j] = 1.0
            elseif Δd <= 0.5
                img[i, j] = max(0.5 - Δd, img[i, j])
            end
        end
    end

    Gray.(img)
end
