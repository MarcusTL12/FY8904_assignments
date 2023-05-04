
function plot_chain(chain, monomer_types)
    xs = [x for (x, _) in chain]
    ys = [y for (_, y) in chain]

    Plots.plot(xs, ys; leg=false, ratio=1.0)
    for (i, monomer) in enumerate(monomer_types)
        Plots.annotate!(xs[i], ys[i], "$monomer")
    end

    Plots.plot!()
end

function plot_chain(chain::AbstractArray{NTuple{2,Int}})
    xs = [x for (x, _) in chain]
    ys = [y for (_, y) in chain]

    Plots.plot(xs, ys; leg=false, ratio=1.0)
    Plots.scatter!(xs, ys)
end

function plot_chain(chain::AbstractArray{NTuple{3,Int}})
    xs = [x for (x, _, _) in chain]
    ys = [y for (_, y, _) in chain]
    zs = [z for (_, _, z) in chain]

    x12, y12, z12 = extrema(xs), extrema(ys), extrema(zs)
    d = max(- -(x12...), - -(y12...), - -(z12...)) / 2
    xm, ym, zm = mean(x12), mean(y12), mean(z12)

    xs = xs .- xm
    ys = ys .- ym
    zs = zs .- zm

    # This was very finicky to get working properly because the aspect
    # ratio options do not work right for 3d. This seems to work nicely.
    Plots.scatter(xs, ys, zs;
        markersize=2, leg=false, size=(600, 600),
        # xlims=(-d, d), ylims=(-d, d), zlims=(-d, d)
    )
    Plots.plot!(xs, ys, zs;
        leg=false, size=(600, 600), linewidth=5,
        # xlims=(-d, d), ylims=(-d, d), zlims=(-d, d)
    )
end

function plot_temperature_interfaces!(p, sweeps)
    Plots.vline!(p, sweeps; color=:black)
end
