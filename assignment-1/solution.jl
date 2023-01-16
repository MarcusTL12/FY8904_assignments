using DataStructures

include("collision-detection.jl")
include("initial-conditions.jl")
include("collision-resolution.jl")
include("visualization.jl")

function test_crater()
    ps, vs, rs, ms = crater_setup(10, 1.0, 1.0, 0.05, 1.0, 1.0)

    pq = init_collisions(ps, vs, rs)

    display(pq)

    (i, j), t1 = dequeue_pair!(pq)

    speed = 300

    img_i = animate_line(ps, vs, rs, 0.0, t1, ceil(Int, speed * t1), 1, "assignment-1/tmp_anim/", 1000)

    ps .+= vs .* t1

    do_collision!(i, j, ps, vs, ms, 1.0)

    update_collisions!(pq, i, t1, ps, vs, rs)
    update_collisions!(pq, j, t1, ps, vs, rs)

    (i, j), t2 = dequeue_pair!(pq)

    animate_line(ps, vs, rs, 0.0, t2 - t1, ceil(Int, speed * (t2 - t1)), img_i, "assignment-1/tmp_anim/", 1000; include_first=false)

    pq
end
