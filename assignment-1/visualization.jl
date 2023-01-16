using Images

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
