using Plots

function calculate_kin_e(vs, ms)
    sum(0.5 * m * v'v for (m, v) in zip(ms, eachrow(vs)))
end

function plot_kin_e_history(t_hist, vs_hist, ms)
    kin_e_hist = [calculate_kin_e(vs, ms) for vs in eachslice(vs_hist; dims=3)]

    rel_kin_e = kin_e_hist ./ kin_e_hist[1]

    plot(t_hist, rel_kin_e; leg=false, ylims=[0.0, 1.0])
end
