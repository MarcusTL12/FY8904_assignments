using DataStructures
using Plots

include("collision-detection.jl")
include("initial-conditions.jl")
include("collision-resolution.jl")
include("visualization.jl")

function test_crater()
    ps, vs, rs, ms = crater_setup(20, 1.0, 1.0, 0.05, 1.0, 1.0)
    # ps, vs, rs, ms = crater_setup(1, 1.0, 1.0, 0.05, 1.0, 1.0)
    # ps, vs, rs, ms = simple()

    ξ = 0.90
    t_target = 5.0

    println("Initializing PQ:")
    pq = @time init_collisions(ps, vs, rs)

    (i, j), t = dequeue_pair!(pq)
    Δt = t

    time_per_frame = 0.003

    anim_dir = "assignment-1/tmp_anim/"
    frames_dir = joinpath(anim_dir, "frames")
    rm(frames_dir; recursive=true, force=true)
    # mkdir(anim_dir)
    mkdir(frames_dir)

    img_i, rem_t = animate_line(
        ps, vs, rs, 0.0, t, time_per_frame, 1,
        frames_dir, 1000
    )

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

        img_i, rem_t = animate_line(
            ps, vs, rs, time_per_frame - rem_t, Δt, time_per_frame, img_i,
            frames_dir, 1000;
            include_first=false
        )

        timer3 = time()

        t_frame_gen += timer3 - timer2

        n_collisions += 1

        if timer3 >= progress_timer + progress_interval
            progress_timer += progress_interval
            println("$(round(t; digits=4)) / $t_target = ",
                round(100 * t / t_target; digits=2), "%")
        end
    end

    make_mp4(frames_dir, anim_dir)

    @show t_logic t_frame_gen n_collisions
end
