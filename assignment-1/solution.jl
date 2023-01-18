using DataStructures
using Plots

include("collision-detection.jl")
include("initial-conditions.jl")
include("collision-resolution.jl")
include("visualization.jl")
include("simulation.jl")
include("analysis.jl")

function test_setup()
    ps, vs, rs, ms = still_hex_grid_box(0.0, 1.0, 0.0, 1.0, 10, 0.5, 1.0)

    show_disks(ps, rs, 1000)
end

function test_crater()
    ps, vs, rs, ms = crater_setup_hex(10, 0.8, 1.0, 0.05, 2.0, 1.0)

    ξ = 0.99
    t_target = 5.0
    kin_e_target = 0.1
    time_per_frame = 0.003
    resolution = 1000

    anim_dir = "assignment-1/tmp_anim/"

    ps_hist, vs_hist, t_hist = simulate(ps, vs, rs, ms, ξ,
        t_target, kin_e_target,
        time_per_frame, true, anim_dir, resolution, true)

    plot_kin_e_history(t_hist, vs_hist, ms)
end

function test_gas()
    ps, vs, rs, ms = hexgrid_rand_velocities(100, 0.01, 1.0, 1.0)

    ξ = 1.0
    t_target = 5.0
    kin_e_target = 0.1
    time_per_frame = 0.003
    resolution = 1000

    anim_dir = "assignment-1/tmp_anim/"

    simulate(ps, vs, rs, ms, ξ,
        t_target, kin_e_target,
        time_per_frame, true, anim_dir, resolution, false)
end
