using DataStructures
using Plots

include("collision-detection.jl")
include("initial-conditions.jl")
include("collision-resolution.jl")
include("visualization.jl")

function test_crater()
    ps, vs, rs, ms = crater_setup(30, 20.0, 1.0, 0.2, 10.0, 1.0)
    # ps, vs, rs, ms = crater_setup(1, 1.0, 1.0, 0.05, 1.0, 1.0)
    # ps, vs, rs, ms = simple()

    pq = init_collisions(ps, vs, rs)

    (i, j), t = dequeue_pair!(pq)
    Δt = t

    time_per_frame = 0.01

    anim_dir = "assignment-1/tmp_anim/"
    frames_dir = joinpath(anim_dir, "frames")
    rm(frames_dir; recursive=true, force=true)
    # mkdir(anim_dir)
    mkdir(frames_dir)

    img_i, rem_t = animate_line(
        ps, vs, rs, 0.0, t, time_per_frame, 1,
        frames_dir, 1000
    )

    while t < 5
        ps .+= vs .* Δt

        do_collision!(i, j, ps, vs, ms, 1.0)

        update_collisions!(pq, i, t, ps, vs, rs)
        update_collisions!(pq, j, t, ps, vs, rs)

        (i, j), nt = dequeue_pair!(pq)

        Δt = nt - t
        t = nt

        img_i, rem_t = animate_line(
            ps, vs, rs, time_per_frame - rem_t, Δt, time_per_frame, img_i,
            frames_dir, 1000;
            include_first=false
        )
    end

    make_mp4(frames_dir, anim_dir)
end
