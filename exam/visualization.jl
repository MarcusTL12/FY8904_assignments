
function plot_chain(chain, monomer_types, annotate=true)
    xs = [x for (x, _) in chain]
    ys = [y for (_, y) in chain]

    Plots.plot(xs, ys; leg=false, ratio=1.0)

    if annotate
        for (i, monomer) in enumerate(monomer_types)
            Plots.annotate!(xs[i], ys[i], "$monomer")
        end
    else
        Plots.scatter!(xs, ys)
    end

    Plots.plot!()
end
