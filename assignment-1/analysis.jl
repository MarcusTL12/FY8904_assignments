using Plots
using Statistics
using KernelDensity
using StatsPlots

function calculate_kin_e(vs, ms)
    sum(0.5 * m * v'v for (m, v) in zip(ms, eachrow(vs)))
end

function plot_kin_e_history(t_hist, vs_hist, ms)
    kin_e_hist = [calculate_kin_e(vs, ms) for vs in eachslice(vs_hist; dims=3)]

    rel_kin_e = kin_e_hist ./ kin_e_hist[1]

    plot(t_hist, rel_kin_e; leg=false, ylims=[0.0, 1.0])
end

function plot_kin_e_history!(t_hist, vs_hist, ms, tot_e; label="")
    kin_e_hist = [calculate_kin_e(vs, ms) for vs in eachslice(vs_hist; dims=3)]

    rel_kin_e = kin_e_hist ./ tot_e

    plot!(t_hist, rel_kin_e; label=label)
end

function calculate_mean_v(vs)
    mean(hypot(v[1], v[2]) for v in eachrow(vs))
end

function calculate_std_v(vs)
    std(hypot(v[1], v[2]) for v in eachrow(vs))
end

function plot_v_mean_stddev(t_hist, vs_hist; label="")
    means = [calculate_mean_v(vs) for vs in eachslice(vs_hist; dims=3)]
    stds = [calculate_std_v(vs) for vs in eachslice(vs_hist; dims=3)]

    plot(t_hist, means; label="mean $label")
    plot!(t_hist, stds; label="std $label")
end

function plot_v_mean_stddev!(t_hist, vs_hist; label="")
    means = [calculate_mean_v(vs) for vs in eachslice(vs_hist; dims=3)]
    stds = [calculate_std_v(vs) for vs in eachslice(vs_hist; dims=3)]

    plot!(t_hist, means; label="mean $label")
    plot!(t_hist, stds; label="std $label")
end

function plot_v_dist(t, t_hist, vs_hist, bandwidth=nothing)
    i = get_closest_t_indices(t_hist, [t])[1]

    vs = @view vs_hist[:,:,i]

    data = [hypot(v[1], v[2]) for v in eachrow(vs)]

    if !isnothing(bandwidth)
        density!(data; label="t=$(round(t; digits=2))", bandwidth=bandwidth, trim=true)
    else
        density!(data; label="t=$(round(t; digits=2))", trim=true)
    end
end

function plot_v_dist_window(tmin, tmax, nt, t_hist, vs_hist, bandwidth=nothing; label="")
    ts = range(tmin, tmax, nt)

    is = get_closest_t_indices(t_hist, ts)

    data = Float64[]

    for i in is
        vs = @view vs_hist[:,:,i]

        append!(data, hypot(v[1], v[2]) for v in eachrow(vs))
    end

    if !isnothing(bandwidth)
        density!(data; label="t∈($(round(tmin; digits=2)),$(round(tmax; digits=2))),n=$nt,$label", bandwidth=bandwidth, trim=true)
    else
        density!(data; label="t∈($(round(tmin; digits=2)),$(round(tmax; digits=2))),n=$nt,$label", trim=true)
    end
end

function plot_v_dist_hist(ts, t_hist, vs_hist)
    is = get_closest_t_indices(t_hist, ts)

    plot()
    for i in is
        t = round(t_hist[i]; digits=2)
        vs = @view vs_hist[:,:,i]

        data = [hypot(v[1], v[2]) for v in eachrow(vs)]

        density!(data; label="t=$t")
    end
    plot!()
end

function get_closest_t_indices(t_hist, ts)
    is = Int[]

    last_i = 1

    for t_target in ts
        did_push = false
        for i in last_i+1:length(t_hist)
            if t_hist[i] > t_target
                j = if abs(t_target - t_hist[i]) < abs(t_target - t_hist[i-1])
                    i
                else
                    i - 1
                end

                push!(is, j)
                last_i = j
                did_push = true
                break
            end
        end

        if !did_push
            push!(is, length(t_hist))
        end
    end

    is
end

function plot_theoretical_maxwell_boltzmann(E, m, n, vmin, vmax; label="")
    kbT = E / n

    vr = range(vmin, vmax, 1000)

    plot!(vr, v -> m * v / kbT * exp(-m * v^2 / 2kbT); label="MBD,$label")
end
