using Images
using SIMD
using LinearAlgebra
using FFMPEG

include("simd-utils.jl")

function show_disks(ps, rs, resolution)
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
                # img[i, j] = max(0.5 - Δd, img[i, j])
                # img[i, j] = 0.5 - Δd
                img[i, j] = min(1.0, img[i, j] + 0.5 - Δd)
            end
        end
    end

    Gray.(img)
end

function get_color(m, m_mid)
    if m < m_mid
        RGB(0.0392156862745098, 0.8392156862745098, 1.0)
    else
        RGB(1.0, 0.0, 0.2)
    end
end

function show_disks_colored(ps, rs, ms, m_mid, resolution)
    img = zeros(RGB{Float64}, resolution, resolution)

    xr = range(0.0, 1.0, resolution)
    yr = range(0.0, 1.0, resolution)

    @inbounds for ((px, py), r, m) in zip(eachrow(ps), rs, ms)
        col = get_color(m, m_mid)

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
                img[i, j] = clamp01(img[i, j] + col)
            elseif Δd <= 0.5
                img[i, j] = clamp01(img[i, j] + (0.5 - Δd) * col)
            end
        end
    end

    img
end

function show_disks_threaded(ps, rs, resolution)
    img = zeros(resolution, resolution)

    xr = range(0.0, 1.0, resolution)
    yr = range(0.0, 1.0, resolution)

    nth = Threads.nthreads()

    Threads.@threads for id in 1:nth
        for k in id:nth:length(rs)
            px = ps[k, 1]
            py = ps[k, 2]
            r = rs[k]

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
    end

    Gray.(img)
end

function animate_line(ps, vs, rs, ms, m_mid, t1, t2, Δt, img_i, dir,
    resolution, color)
    tr = t1:Δt:t2

    for t in tr
        if color
            img = show_disks_colored(ps, rs, ms, m_mid, resolution)
        else
            img = show_disks(ps .+ vs .* t, rs, resolution)
        end

        save("$dir/$(lpad(img_i, 6, "0")).png", img)

        img_i += 1
    end

    img_i, t2 - last(tr)
end

function make_mp4(framesdir, outdir; nth=Threads.nthreads())
    FFMPEG.exe(`-threads $nth -framerate 60 -i $framesdir/%06d.png -c:v libx264 -pix_fmt yuv420p $outdir/out.mp4 -y`)
end

function compress_output(outdir; nth=Threads.nthreads())
    FFMPEG.exe(`-threads $nth -i $outdir/out.mp4 -vcodec libx265 -crf 28 $outdir/out_cmp.mp4`)
end
