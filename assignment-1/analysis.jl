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

function calculate_mean_v(vs)
    mean(hypot(v[1], v[2]) for v in eachrow(vs))
end

function calculate_std_v(vs)
    std(hypot(v[1], v[2]) for v in eachrow(vs))
end

function plot_v_mean_stddev(t_hist, vs_hist)
    means = [calculate_mean_v(vs) for vs in eachslice(vs_hist; dims=3)]
    stds = [calculate_std_v(vs) for vs in eachslice(vs_hist; dims=3)]

    plot(t_hist, means; label="mean")
    plot!(t_hist, stds; label="std")
end

function plot_v_dist(vs)
    data = [hypot(v[1], v[2]) for v in eachrow(vs)]

    k = kde(data)

    plot(0:0.01:5, x -> pdf(k, x); leg=false)
end

function plot_v_dist_hist(ts, t_hist, vs_hist, max_v)
    is = get_closest_t_indices(t_hist, ts)

    plot()
    for i in is
        t = round(t_hist[i]; digits=2)
        vs = @view vs_hist[:,:,i]

        data = [hypot(v[1], v[2]) for v in eachrow(vs)]

        # k = kde(data)

        # plot!(range(0.0, max_v, 1000), x -> pdf(k, x); label="t=$t")

        # histogram!(data; label="t=$t", bins=100)
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
