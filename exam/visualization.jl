
# This function plots a 2d polymer chain and annotates all the monomers by
# the amino acid type (number 1 to 20)
function plot_chain(chain, monomer_types)
    xs = [x for (x, _) in chain]
    ys = [y for (_, y) in chain]

    Plots.plot(xs, ys; leg=false, ratio=1.0)
    for (i, monomer) in enumerate(monomer_types)
        Plots.annotate!(xs[i], ys[i], "$monomer")
    end

    Plots.plot!()
end

# This function plots a 2d polymer chain without annotating the monomer types
function plot_chain(chain::AbstractArray{NTuple{2,Int}})
    xs = [x for (x, _) in chain]
    ys = [y for (_, y) in chain]

    Plots.plot(xs, ys; leg=false, ratio=1.0)
    Plots.scatter!(xs, ys)
end

# This plots a 3d chain
function plot_chain(chain::AbstractArray{NTuple{3,Int}})
    xs = [x for (x, _, _) in chain]
    ys = [y for (_, y, _) in chain]
    zs = [z for (_, _, z) in chain]

    x12, y12, z12 = extrema(xs), extrema(ys), extrema(zs)
    xm, ym, zm = mean(x12), mean(y12), mean(z12)

    xs = xs .- xm
    ys = ys .- ym
    zs = zs .- zm

    # This was very finicky to get working properly because the aspect
    # ratio options do not work right for 3d.
    # There is still an open issue about this from 2019:
    # https://github.com/JuliaPlots/Plots.jl/issues/1949
    # This seems to work most of the time.
    Plots.scatter(xs, ys, zs;
        markersize=2, leg=false, size=(600, 600),
    )
    Plots.plot!(xs, ys, zs;
        leg=false, size=(600, 600), linewidth=5,
    )
end

function plot_temperature_interfaces!(p, sweeps)
    Plots.vline!(p, sweeps; color=:black)
end
