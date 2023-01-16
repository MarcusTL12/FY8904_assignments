using DataStructures

include("collision-detection.jl")
include("initial-conditions.jl")
include("collision-resolution.jl")
include("visualization.jl")

function test_crater()
    ps, vs, rs, ms = crater_setup(10, 1.0, 1.0, 0.05, 1.0, 1.0)

    pq = init_collisions(ps, vs, rs)

    (i, j), t = dequeue_pair!(pq)

    ps .+= vs .* t

    show_disks(ps, rs, 1000)
end
