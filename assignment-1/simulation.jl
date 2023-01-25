
function simulate(ps, vs, rs, ms, ξ, t_target, kin_e_target, collission_target, time_per_frame,
    animate, color, anim_dir, resolution, collect_data, data_interval)
    println("Initializing PQ:")
    pq = @time init_collisions(ps, vs, rs)

    m_mid = (ms[1] + ms[end]) / 2

    ps_history = copy(ps[:])
    vs_history = copy(vs[:])
    t_history = [0.0]

    init_kin_e = calculate_kin_e(vs, ms)

    (i, j), t = dequeue_pair!(pq)
    Δt = t

    if animate
        frames_dir = joinpath(anim_dir, "frames")
        rm(frames_dir; recursive=true, force=true)
        mkdir(frames_dir)

        img_i, rem_t = animate_line(
            ps, vs, rs, ms, m_mid, 0.0, t, time_per_frame, 1,
            frames_dir, resolution, color
        )
    end

    t_logic = 0.0
    t_frame_gen = 0.0

    n_collisions = 1
    n_data_points = 1

    progress_timer = time()
    progress_interval = 2.0

    while t < t_target &&
              calculate_kin_e(vs, ms) / init_kin_e > kin_e_target &&
              n_collisions < collission_target
        timer1 = time()
        ps .+= vs .* Δt

        do_collision!(i, j, ps, vs, ms, ξ)

        update_collisions!(pq, i, t, ps, vs, rs)
        update_collisions!(pq, j, t, ps, vs, rs)

        if collect_data && n_collisions % data_interval == 0
            append!(ps_history, ps)
            append!(vs_history, vs)
            push!(t_history, t)
            n_data_points += 1
        end

        (i, j), nt = dequeue_pair!(pq)

        Δt = nt - t
        t = nt

        timer2 = time()

        t_logic += timer2 - timer1

        if animate
            img_i, rem_t = animate_line(
                ps, vs, rs, ms, m_mid, time_per_frame - rem_t, Δt,
                time_per_frame, img_i, frames_dir, resolution, color
            )
        end

        timer3 = time()

        t_frame_gen += timer3 - timer2

        n_collisions += 1

        if timer3 >= progress_timer + progress_interval
            progress_timer += progress_interval
            println("$n_collisions / $collission_target = ",
                round(100 * n_collisions / collission_target; digits=2), "%")
            println("$(round(t; digits=4)) / $t_target = ",
                round(100 * t / t_target; digits=2), "%")
            println("Ek/Ek0: ",
                round(100 * calculate_kin_e(vs, ms) / init_kin_e; digits=2),
                "% > ", round(100 * kin_e_target; digits=2), "%\n"
            )
        end
    end

    if collect_data
        append!(ps_history, ps)
        append!(vs_history, vs)
        push!(t_history, t)
        n_data_points += 1
    end

    if animate
        make_mp4(frames_dir, anim_dir)
    end

    @show t_logic t_frame_gen n_collisions

    if collect_data
        reshape(ps_history, size(ps)..., n_data_points),
        reshape(vs_history, size(vs)..., n_data_points),
        t_history
    end
end
