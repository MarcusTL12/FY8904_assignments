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
    ps, vs, rs, ms = hexgrid_rand_velocities(50, 0.01, 1.0, 1.0)

    E = calculate_kin_e(vs, ms)
    n = size(ps, 1)
    m = ms[1]

    @show n

    ξ = 1.0
    t_target = 10.0
    kin_e_target = 0.1
    time_per_frame = 0.003
    resolution = 1000

    anim_dir = "assignment-1/tmp_anim/"

    ps_hist, vs_hist, t_hist = simulate(ps, vs, rs, ms, ξ,
        t_target, kin_e_target,
        time_per_frame, false, anim_dir, resolution, true)

    display(plot_v_mean_stddev(t_hist, vs_hist))
    plot()
    plot_v_dist_window(0.5, t_target, 1000, t_hist, vs_hist)
    plot_theoretical_maxwell_boltzmann(E, m, n, 0.0, 3.0)
end

function test_inhomogenous_gas()
    ps, vs, rs, ms = hexgrid_rand_velocities(40, 0.1, 1.0, 1.0)

    n = size(ps, 1)

    @show n

    n_half = n ÷ 2
    is_lo = 1:n_half
    is_hi = n_half+1:n

    ms[n_half + 1:end] *= 4

    E = calculate_kin_e(vs, ms)
    m_lo = ms[1]
    m_hi = ms[end]

    ξ = 1.0
    t_target = 10.0
    kin_e_target = 0.1
    time_per_frame = 0.003
    resolution = 1000

    anim_dir = "assignment-1/tmp_anim/"

    ps_hist, vs_hist, t_hist = simulate(ps, vs, rs, ms, ξ,
        t_target, kin_e_target,
        time_per_frame, false, anim_dir, resolution, true)

    t_start = 5.0
    n_points = 1000

    plot_v_mean_stddev(t_hist, vs_hist[is_lo,:,:]; label="m=1")
    display(plot_v_mean_stddev!(t_hist, vs_hist[is_hi,:,:]; label="m=4"))
    plot()
    plot_v_dist_window(t_start, t_target, n_points, t_hist, vs_hist[is_lo,:,:]; label="m=1")
    plot_theoretical_maxwell_boltzmann(E / 2, m_lo, n_half, 0.0, 5.0; label="m=1")
    plot_v_dist_window(t_start, t_target, n_points, t_hist, vs_hist[is_hi,:,:]; label="m=4")
    plot_theoretical_maxwell_boltzmann(E / 2, m_hi, n_half, 0.0, 3.0; label="m=4")
end
