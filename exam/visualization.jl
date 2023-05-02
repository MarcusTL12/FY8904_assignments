
function plot_chain(chain, monomer_types)
    xs = [x for (x, _) in chain]
    ys = [y for (_, y) in chain]

    Plots.plot(xs, ys; leg=false, ratio=1.0)

    for (i, monomer) in enumerate(monomer_types)
        Plots.annotate!(xs[i], ys[i], "$monomer")
    end

    Plots.plot!()
end
