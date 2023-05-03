
function plot_chain(chain, monomer_types)
    xs = [x for (x, _) in chain]
    ys = [y for (_, y) in chain]

    Plots.plot(xs, ys; leg=false, ratio=1.0)
    for (i, monomer) in enumerate(monomer_types)
        Plots.annotate!(xs[i], ys[i], "$monomer")
    end

    Plots.plot!()
end

function plot_chain(chain)
    xs = [x for (x, _) in chain]
    ys = [y for (_, y) in chain]

    Plots.plot(xs, ys; leg=false, ratio=1.0)
    Plots.scatter!(xs, ys)

    Plots.plot!()
end

function plot_temperature_interfaces!(p, sweeps)
    Plots.vline!(p, sweeps; color=:black)
end
