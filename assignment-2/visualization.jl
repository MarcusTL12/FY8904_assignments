
function animate_spin_history(lattice_points, spin_history)
    f = Figure(resolution=(3840, 2160))
    ax = Axis3(f[1, 1], viewmode=:fitzoom, aspect=:data, perspectiveness=0.5)

    spinx = Observable(@view spin_history[:, 1, 1])
    spiny = Observable(@view spin_history[:, 2, 1])
    spinz = Observable(@view spin_history[:, 3, 1])

    arrows!(ax,
        (@view lattice_points[:, 1]),
        (@view lattice_points[:, 2]),
        (@view lattice_points[:, 3]),
        spinx, spiny, spinz;
        lengthscale=1.0, linewidth=0.05, arrowsize=0.1, color=:gray
    )

    record(f, "tmp.mp4", axes(spin_history, 3), framerate=60) do i
        spinx[] = @view spin_history[:, 1, i]
        spiny[] = @view spin_history[:, 2, i]
        spinz[] = @view spin_history[:, 3, i]
    end
end

function visualize_spin_history_interactive(lattice_points, spin_history)
    f = Figure(resolution=(1500, 900))
    ax = Axis3(f[1, 1], viewmode=:fitzoom, aspect=:data, perspectiveness=0.5)

    xmin, xmax = extrema(@view lattice_points[:, 1])
    ymin, ymax = extrema(@view lattice_points[:, 2])
    zmin, zmax = extrema(@view lattice_points[:, 3])

    limits!(ax,
        xmin - 1, xmax + 1,
        ymin - 1, ymax + 1,
        zmin - 1, zmax + 1
    )

    framerate = 60.0
    frametime = 1.0 / framerate

    frame_n = Observable(1)

    spinx = @lift (@view spin_history[:, 1, $frame_n])
    spiny = @lift (@view spin_history[:, 2, $frame_n])
    spinz = @lift (@view spin_history[:, 3, $frame_n])

    arrows!(ax,
        (@view lattice_points[:, 1]),
        (@view lattice_points[:, 2]),
        (@view lattice_points[:, 3]),
        spinx, spiny, spinz;
        lengthscale=0.4, linewidth=0.04, arrowsize=0.1, color=:gray,
        arrowcolor=:red
    )

    isrunning = Observable(false)

    sl = Slider(f[2, 1][1, 1], range=axes(spin_history, 3))

    lift(sl.value) do i
        frame_n[] = i
    end

    play_label = @lift $isrunning ? "▮▮" : " ▸ "
    play = Button(f[2, 1][1, 2]; label=play_label)

    on(play.clicks) do _
        if !isrunning[]
            isrunning[] = true
            @async begin
                t = time()
                if frame_n[] == size(spin_history, 3)
                    frame_n[] = 1
                end

                while frame_n[] < size(spin_history, 3)
                    if !isrunning[]
                        break
                    end

                    target_time = t + frametime
                    curtime = time()
                    sleep(max(target_time - curtime, 0.0))
                    t = target_time

                    frame_n[] += 1
                    set_close_to!(sl, frame_n[])
                end

                isrunning[] = false
            end
        else
            isrunning[] = false
        end
    end

    reset = Button(f[2, 1][1, 3]; label="⬛")

    on(reset.clicks) do _
        @async begin
            isrunning[] = false
            frame_n[] = 1
            set_close_to!(sl, frame_n[])
        end
    end

    rec_btn = Button(f[2, 1][1, 4]; label="●")

    on(rec_btn.clicks) do _
        @async record(f, "tmp.mp4", axes(spin_history, 3), framerate=60) do i
            frame_n[] = i
            set_close_to!(sl, frame_n[])
        end
    end

    f
end

function spin_to_HSL(sx, sy, sz)
    HSL(atand(sx, sy), 1.0, (sz + 1.0) * 0.5)
end

function visualize_spin_surface!(S)
    @inbounds [RGB(spin_to_HSL(S[x, y, 1], S[x, y, 2], S[x, y, 3]))
               for x in axes(S, 1), y in axes(S, 1)]
end

function animate_2d_spin_history(S_hist, nx, ny)
    S_hist = reshape(S_hist, nx, ny, 3, size(S_hist, 3))

    dir = "tmp_frames"

    rm(dir; recursive=true, force=true)
    mkdir(dir)

    @time for i in axes(S_hist, 4)
        img = visualize_spin_surface(@view S_hist[:, :, :, i])

        save("$dir/$i.png", img)
    end

    make_mp4(dir, "tmp_anim")
end

function make_mp4(framesdir, outdir; nth=Threads.nthreads())
    FFMPEG.exe(`-threads $nth -framerate 60 -i $framesdir/%d.png -c:v libx264 \
    -pix_fmt yuv420p $outdir/out.mp4 -y`)
end
