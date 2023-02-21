using GLMakie

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
    f = Figure(resolution=(1100, 1000))
    ax = Axis3(f[1, 1], viewmode=:fitzoom, aspect=:data, perspectiveness=0.5)

    frame_n = Observable(1)

    spinx = @lift (@view spin_history[:, 1, $frame_n])
    spiny = @lift (@view spin_history[:, 2, $frame_n])
    spinz = @lift (@view spin_history[:, 3, $frame_n])

    arrows!(ax,
        (@view lattice_points[:, 1]),
        (@view lattice_points[:, 2]),
        (@view lattice_points[:, 3]),
        spinx, spiny, spinz;
        lengthscale=0.4, linewidth=0.04, arrowsize=0.1, color=:gray
    )

    isrunning = Observable(false)

    sl = Slider(f[2, 1], range=axes(spin_history, 3))

    lift(sl.value) do i
        @async frame_n[] = i
    end

    play_label = @lift $isrunning ? "▋▋" : "▶"
    play = Button(f[3, 2]; label=play_label)

    on(play.clicks) do _
        if !isrunning[]
            isrunning[] = true
            @async begin
                if frame_n[] == size(spin_history, 3)
                    frame_n[] = 1
                end

                while frame_n[] < size(spin_history, 3)
                    if !isrunning[]
                        break
                    end
                    frame_n[] += 1
                    set_close_to!(sl, frame_n[])
                    sleep(1 / 60)
                end
                isrunning[] = false
            end
        else
            isrunning[] = false
        end
    end

    reset = Button(f[4, 2]; label="█▋")

    on(reset.clicks) do _
        isrunning[] = false
        frame_n[] = 1
    end

    f
end
