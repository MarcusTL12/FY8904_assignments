using DataStructures
using Plots

include("collision-detection.jl")
include("initial-conditions.jl")
include("collision-resolution.jl")
include("visualization.jl")
include("simulation.jl")

function test_crater()
    ps, vs, rs, ms = crater_setup(50, 0.1, 1.0, 0.05, 1.0, 1.0)

    ξ = 1.0
    t_target = 5.0
    time_per_frame = 0.003
    resolution = 1000

    anim_dir = "assignment-1/tmp_anim/"

    simulate(ps, vs, rs, ms, ξ, t_target, time_per_frame, true, anim_dir, resolution)
end
