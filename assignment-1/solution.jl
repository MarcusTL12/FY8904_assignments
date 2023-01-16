using DataStructures

include("collision-detection.jl")
include("initial-conditions.jl")
include("collision-resolution.jl")
include("visualization.jl")

function test_crater()
    ps, vs, rs, ms = crater_setup(10, 1.0, 1.0, 0.05, 1.0, 1.0)

    pq = init_collisions(ps, vs, rs)

    (i, j), t = dequeue_pair!(pq)

    img_i = animate_line(ps, vs, rs, 0.0, t, 100, 1, "assignment-1/tmp_anim/", 1000)

    ps .+= vs .* t

    do_collision!(i, j, ps, vs, ms, 1.0)

    animate_line(ps, vs, rs, 0.0, 0.05, 30, img_i, "assignment-1/tmp_anim/", 1000; include_first=false)
end
