using LinearAlgebra
import Plots
using Printf
using DSP

include("visualization.jl")
include("lattice2d.jl")
include("mc2d.jl")

# This is run to "solve" task 2.1.3
# Here we test a few different foldings of a 15 monomer long protein
function run_2_1_3()
    # filepath to save the figures
    figure_path = "2.1.3-figures"

    # Initialize random interaction energies
    interaction_matrix = make_interaction_energy_matrix()
    p = Plots.heatmap(Matrix(interaction_matrix)[:, end:-1:1];
        ratio=1.0, margins=0.2 * Plots.cm)
    display(p)
    Plots.savefig(p, joinpath(figure_path, "matrix.pdf"))

    # Choosing a polymer length and initializing the monomer types as random
    # integers in the range 1:20
    n = 15
    monomer_types = rand(1:20, n)

    # Here we test different foldings. The 'width' w indicate how many
    # monmers the sheet is wide. This means that both 15 and 1 should be
    # fully unfolded, just along different axes. w = 5 or 3 are special as the
    # 15 length polymer fits in a perfect 3x5 or 5x3 rectangle, which maximizes
    # the number of interactions.
    for w in (10, 7, 5, 3, 2, 1)
        chain = make_folded_2d_chain(n, w)
        coord_map = make_coord_map(chain)
        energy = calculate_energy_direct(interaction_matrix, monomer_types,
            chain, coord_map)

        n_interactions = length(make_neighbour_list(chain, coord_map))

        p = plot_chain(chain, monomer_types)
        Plots.savefig(p, joinpath(figure_path, "$w.pdf"))

        @printf "w = %d: %.2f\n" w energy
        @printf "Number of interactions: %d\n\n" n_interactions
    end
end

function run_2_1_5()
    figure_path = "2.1.5-figures"

    interaction_matrix = make_interaction_energy_matrix()

    n = 15
    monomer_types = rand(1:20, n)

    chain = make_linear_2d_chain(n)
    coord_map = make_coord_map(chain)

    snapshot_inds = (1, 10, 100)

    energies, e2e_dists, RoGs, snapshots = @time simulate_2d(
        chain, coord_map, interaction_matrix, monomer_types, 10, 100,
        snapshot_inds
    )

    Plots.savefig(Plots.plot(energies), joinpath(figure_path, "energy.pdf"))
    Plots.savefig(Plots.plot(e2e_dists), joinpath(figure_path, "e2e.pdf"))
    Plots.savefig(Plots.plot(RoGs), joinpath(figure_path, "RoG.pdf"))

    for (snap, i) in zip(eachcol(snapshots), snapshot_inds)
        Plots.savefig(plot_chain(snap), joinpath(figure_path, "snapshot$i.pdf"))
    end
end

function calc_window_avg(xs, n)
    w = [1 / n for _ in 1:n]

    conv(xs, w)[1:length(xs)]
end

function test_mc2d()
    interaction_matrix = make_interaction_energy_matrix()

    n = 100
    monomer_types = rand(1:20, n)

    chain = make_linear_2d_chain(n)
    coord_map = make_coord_map(chain)

    energies, e2e_dists, RoGs = @time simulate_2d(
        chain, coord_map, interaction_matrix, monomer_types, 5, 100000
    )

    display(Plots.plot(energies; title="Energy"))
    display(Plots.plot(calc_window_avg(energies, 1000); title="Energy avg"))
    display(Plots.plot(e2e_dists; title="End to end dist"))
    display(Plots.plot(RoGs; title="Radius of gyration"))

    plot_chain(chain)
end
