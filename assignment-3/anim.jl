using GLMakie
using Printf
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

    gain_sl = Slider(f[1, 2][1, 1][1, 1], range=0.5:0.01:2.0, startvalue=1.0,
        horizontal=false).value
    Label(f[1, 2][1, 1][2, 1],
        @lift @sprintf "Gain: ↑ = %.2f\nFlip sign: ↓" $gain_sl)
    sign_tgl = Toggle(f[1, 2][1, 1][3, 1]).active
    gain = @lift $sign_tgl ? -$gain_sl : $gain_sl
    scaled_mode = Observable(copy(unpacked_mode[]))

    speed_sl = Slider(f[1, 2][1, 1][1, 2], range=0.0:0.01:10.0, startvalue=1.0,
        horizontal=false).value
    Label(f[1, 2][1, 1][2, 2], @lift @sprintf "Speed: %.2f" $speed_sl)

    on(unpacked_mode) do mode
        smode = scaled_mode[]
        g = gain[]
        for (i, x) in enumerate(mode)
            smode[i] = x * g
        end
        scaled_mode[] = smode
    end

    on(gain) do g
        smode = scaled_mode[]
        for (i, x) in enumerate(unpacked_mode[])
            smode[i] = x * g
        end
        scaled_mode[] = smode
    end

    t = Observable(0.0)
    time_evolving_mode = Observable(copy(scaled_mode[]))

    on(t) do t
        tmode = time_evolving_mode[]
        ω = √(e[mode_i[]])
        for (i, x) in enumerate(scaled_mode[])
            tmode[i] = x * cos(ω * t)
        end
        time_evolving_mode[] = tmode
    end

    on(scaled_mode) do m
        tmode = time_evolving_mode[]
        ω = √(e[mode_i[]])
        t_ = t[]
        for (i, x) in enumerate(m)
            tmode[i] = x * cos(ω * t_)
        end
        time_evolving_mode[] = tmode
    end

    GLMakie.surface!(ax3d, xy_range, xy_range, time_evolving_mode;
        colorrange=(@lift $colorrange .* $gain_sl))
    hm = GLMakie.heatmap!(ax2d, xy_range, xy_range, unpacked_mode;
        colorrange=colorrange)

    Colorbar(f[1, 1][2, 1][1, 2], hm)

    on(mode_i) do i
        mode = unpacked_mode[]
        unpack_mode!(mode, lattice, (@view v[:, i]))
        unpacked_mode[] = mode
    end

    inc_button = Button(config_panel[1, 1][1, 1], label="▲")
    dec_button = Button(config_panel[1, 1][1, 2], label="▼")
    _ = Label(config_panel[1, 2],
        (@lift begin
            @sprintf "Mode: %i/%i \n\
ω/v = %.1f" $mode_i length(e) √(e[$mode_i])
        end))
    mode_box = Textbox(config_panel[2, 2], validator=Int,
        placeholder=":")
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

    isrunning = Observable(false)
    play_label = @lift $isrunning ? "▮▮" : " ▸ "
    play_button = Button(config_panel[2, 1][1, 1], label=play_label)
    reset_button = Button(config_panel[2, 1][1, 2]; label="⬛")

    playback_speed = @lift $speed_sl * 2π / e[1]
    framerate = 60.0
    frametime = 1.0 / framerate

    on(play_button.clicks) do _
        if !isrunning[]
            isrunning[] = true
            @async begin
                timer = time()

                while isrunning[]
                    target_time = timer + frametime
                    curtime = time()
                    sleep(max(target_time - curtime, 0.0))
                    timer = target_time

                    t[] += playback_speed[] * frametime
                end
            end
        else
            isrunning[] = false
        end
    end

    on(reset_button.clicks) do _
        isrunning[] = false
        t[] = 0.0
    end

    f
end
