using GLMakie
GLMakie.activate!(title="Assignment-3", framerate=60.0)

function animate_modes(e, v, lattice, h)
    f = Figure(resolution=(1100, 1600))

    max_z_displacement = maximum(abs, v)

    xy_range = range(0.0; length=size(lattice, 1), step=h)

    aspect = (1.0, 1.0, 0.5)

    config_panel = f[1, 2][2, 1]

    ax3d = Axis3(f[1, 1][1, 1],
        viewmode=:fitzoom, aspect=aspect, perspectiveness=0.5)
    ax2d = Axis(f[1, 1][2, 1][1, 1], aspect=1)

    limits!(ax3d,
        0.0, xy_range[end],
        0.0, xy_range[end],
        -max_z_displacement, max_z_displacement,
    )

    mode_i = Observable(1)

    cur_max_z_displacement = @lift maximum(abs, (@view v[:, $mode_i]))
    colorrange = @lift (-$cur_max_z_displacement, $cur_max_z_displacement)

    unpacked_mode = Observable(unpack_mode(lattice, (@view v[:, 1])))
    GLMakie.surface!(ax3d, xy_range, xy_range, unpacked_mode;
        colorrange=colorrange)
    hm = GLMakie.heatmap!(ax2d, xy_range, xy_range, unpacked_mode;
        colorrange=colorrange)
    
    Colorbar(f[1, 1][2, 1][1, 2], hm)

    on(mode_i) do i
        mode = unpacked_mode[]
        unpack_mode!(mode, lattice, (@view v[:, i]))
        unpacked_mode[] = mode
    end

    inc_button = Button(config_panel[1, 1], label="▲")
    dec_button = Button(config_panel[2, 1], label="▼")
    _ = Label(config_panel[1, 2], (@lift "Mode: $($mode_i)"))
    mode_box = Textbox(config_panel[2, 2], validator=Int,
        placeholder="Enter mode number")
    mode_box.stored_string[] = "1"

    on(mode_box.stored_string) do s
        mode_i[] = parse(Int, s)
    end

    on(inc_button.clicks) do _
        if mode_i[] < length(e)
            mode_i[] += 1
        end
    end

    on(dec_button.clicks) do _
        if mode_i[] > 1
            mode_i[] -= 1
        end
    end

    f
end
