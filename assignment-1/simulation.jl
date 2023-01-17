
function simulate(ps, vs, rs, ms, ξ, t_target, time_per_frame, animate, anim_dir, resolution)
    println("Initializing PQ:")
    pq = @time init_collisions(ps, vs, rs)

    (i, j), t = dequeue_pair!(pq)
    Δt = t

    if animate
        frames_dir = joinpath(anim_dir, "frames")
        rm(frames_dir; recursive=true, force=true)
        mkdir(frames_dir)

        img_i, rem_t = animate_line(
            ps, vs, rs, 0.0, t, time_per_frame, 1,
            frames_dir, resolution
        )
    end

    t_logic = 0.0
    t_frame_gen = 0.0

    n_collisions = 0

    progress_timer = time()
    progress_interval = 2.0

    while t < t_target
        timer1 = time()
        ps .+= vs .* Δt

        do_collision!(i, j, ps, vs, ms, ξ)

        update_collisions!(pq, i, t, ps, vs, rs)
        update_collisions!(pq, j, t, ps, vs, rs)

        (i, j), nt = dequeue_pair!(pq)

        Δt = nt - t
        t = nt

        timer2 = time()

        t_logic += timer2 - timer1

        if animate
            img_i, rem_t = animate_line(
                ps, vs, rs, time_per_frame - rem_t, Δt, time_per_frame, img_i,
                frames_dir, resolution;
                include_first=false
            )
        end

        timer3 = time()

        t_frame_gen += timer3 - timer2

        n_collisions += 1

        if timer3 >= progress_timer + progress_interval
            progress_timer += progress_interval
            println("$(round(t; digits=4)) / $t_target = ",
                round(100 * t / t_target; digits=2), "%")
        end
    end

    if animate
        make_mp4(frames_dir, anim_dir)
    end

    @show t_logic t_frame_gen n_collisions
end
