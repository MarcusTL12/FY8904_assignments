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
        t_target, kin_e_target, Inf,
        time_per_frame, true, false, anim_dir, resolution, true, 10)

    plot_kin_e_history(t_hist, vs_hist, ms)
end

function test_gas()
    ps, vs, rs, ms = hexgrid_rand_velocities(43, 0.5, 1.0, 1.0)

    E = calculate_kin_e(vs, ms)
    n = size(ps, 1)
    m = ms[1]

    @show n

    ξ = 1.0
    t_target = 10.0
    kin_e_target = 0.1
    collision_target = 100n
    time_per_frame = 0.003
    resolution = 1000

    data_target = 10000
    data_interval = cld(collision_target, data_target)
    @show data_interval

    anim_dir = "assignment-1/tmp_anim/"

    ps_hist, vs_hist, t_hist = simulate(ps, vs, rs, ms, ξ,
        t_target, kin_e_target, collision_target,
        time_per_frame, false, false, anim_dir, resolution, true, data_interval)

    t_end = t_hist[end]
    t_start = t_end / 2
    n_points = 1000

    display(plot_v_mean_stddev(t_hist, vs_hist))
    plot()
    plot_v_dist_window(t_start, t_end, n_points, t_hist, vs_hist)
    plot_theoretical_maxwell_boltzmann(E, m, n, 0.0, 3.0)
end

function test_inhomogenous_gas()
    ps, vs, rs, ms = hexgrid_rand_velocities(50, 0.5, 1.0, 1.0)

    n = size(ps, 1)

    @show n

    n_half = n ÷ 2
    is_lo = 1:n_half
    is_hi = n_half+1:n

    ms[n_half+1:end] *= 4

    E = calculate_kin_e(vs, ms)
    m_lo = ms[1]
    m_hi = ms[end]

    ξ = 1.0
    t_target = Inf
    kin_e_target = 0.1
    collision_target = 100n
    time_per_frame = 0.003
    resolution = 1000

    data_target = 10000
    data_interval = cld(collision_target, data_target)
    @show data_interval

    anim_dir = "assignment-1/tmp_anim/"

    ps_hist, vs_hist, t_hist = simulate(ps, vs, rs, ms, ξ,
        t_target, kin_e_target, collision_target,
        time_per_frame, true, true, anim_dir, resolution, true, data_interval)

    plot()
    plot_kin_e_history!(t_hist, (@view vs_hist[is_lo, :, :]),
        (@view ms[is_lo]), E; label="m = 1")
    plot_kin_e_history!(t_hist, (@view vs_hist[is_hi, :, :]),
        (@view ms[is_hi]), E; label="m = 4")

    display(plot!())

    t_end = t_hist[end]
    t_start = t_end / 2
    n_points = 1000

    plot()
    plot_v_dist_window(t_start, t_end, n_points, t_hist,
        (@view vs_hist[is_lo, :, :]); label="m=1")
    plot_theoretical_maxwell_boltzmann(E / 2, m_lo, n_half, 0.0, 5.0;
        label="m=1")
    plot_v_dist_window(t_start, t_end, n_points, t_hist,
        (@view vs_hist[is_hi, :, :]); label="m=4")
    plot_theoretical_maxwell_boltzmann(E / 2, m_hi, n_half, 0.0, 3.0;
        label="m=4")
end

function test_inhomogenous_gas_inelastic()
    ps, vs, rs, ms = hexgrid_rand_velocities(20, 0.1, 1.0, 1.0)

    n = size(ps, 1)

    @show n

    n_half = n ÷ 2
    is_lo = 1:n_half
    is_hi = n_half+1:n

    ms[n_half+1:end] *= 4

    E = calculate_kin_e(vs, ms)
    m_lo = ms[1]
    m_hi = ms[end]

    ξ = 1.01
    t_target = Inf
    kin_e_target = 0.1
    collision_target = 100n
    time_per_frame = 0.003
    resolution = 1000

    data_target = 1000
    data_interval = cld(collision_target, data_target)
    @show data_interval

    anim_dir = "assignment-1/tmp_anim/"

    ps_hist, vs_hist, t_hist = simulate(ps, vs, rs, ms, ξ,
        t_target, kin_e_target, collision_target,
        time_per_frame, true, true, anim_dir, resolution, true, data_interval)


    plot()
    kin_e_hist = [calculate_kin_e(vs, ms) for vs in eachslice(vs_hist; dims=3)]
    plot_kin_e_history!(t_hist, (@view vs_hist[is_lo, :, :]),
        (@view ms[is_lo]), kin_e_hist; label="m = 1")
    plot_kin_e_history!(t_hist, (@view vs_hist[is_hi, :, :]),
        (@view ms[is_hi]), kin_e_hist; label="m = 4")
    display(plot!())

    plot()
    plot_kin_e_history!(t_hist, vs_hist, ms, E; label="tot")
    plot_kin_e_history!(t_hist, (@view vs_hist[is_lo, :, :]),
        (@view ms[is_lo]), E; label="m = 1")
    plot_kin_e_history!(t_hist, (@view vs_hist[is_hi, :, :]),
        (@view ms[is_hi]), E; label="m = 4")
end
