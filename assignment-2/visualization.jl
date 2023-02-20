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
