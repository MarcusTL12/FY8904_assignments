function setup_distributed(n_workers)
    w = workers()
    if length(w) < n_workers
        addprocs(n_workers - length(w))
        @everywhere include("solution.jl")
    elseif length(w) > n_workers
        Δ = length(w) - n_workers
        rmprocs(w[end-(Δ-1):end])
    end
end

# This is run to "solve" task 2.1.3
# Here we test a few different foldings of a 15 monomer long protein
function run_2_1_3()
    # filepath to save the figures
    figure_path = "figures/2.1.3"

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
    for w in (1, 2, 3, 5, 7, 10)
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
    figure_path = "figures/2.1.5"

    interaction_matrix = make_interaction_energy_matrix()

    n = 15
    monomer_types = rand(1:20, n)

    chain = make_linear_2d_chain(n)
    coord_map = make_coord_map(chain)

    snapshot_inds = (1, 10, 100)

    energies, e2e_dists, RoGs, snapshots = @time simulate(
        chain, coord_map, interaction_matrix, monomer_types, 10, 100,
        snapshot_inds
    )

    Plots.savefig(Plots.plot(energies; leg=false),
        joinpath(figure_path, "energy.pdf"))
    Plots.savefig(Plots.plot(e2e_dists; leg=false),
        joinpath(figure_path, "e2e.pdf"))
    Plots.savefig(Plots.plot(RoGs; leg=false),
        joinpath(figure_path, "RoG.pdf"))

    for (snap, i) in zip(eachcol(snapshots), snapshot_inds)
        Plots.savefig(plot_chain(snap), joinpath(figure_path, "snapshot$i.pdf"))
    end
end

function run_2_1_6()
    figure_path = "figures/2.1.6"

    interaction_matrix = make_interaction_energy_matrix()

    n = 15
    monomer_types = rand(1:20, n)

    chain = make_linear_2d_chain(n)
    coord_map = make_coord_map(chain)

    snapshot_inds = (100, 1000, 3000)

    energies, e2e_dists, RoGs, snapshots = @time simulate(
        chain, coord_map, interaction_matrix, monomer_types, 1, 3000,
        snapshot_inds
    )

    Plots.savefig(Plots.plot(energies; leg=false),
        joinpath(figure_path, "energy.pdf"))
    Plots.savefig(Plots.plot(e2e_dists; leg=false),
        joinpath(figure_path, "e2e.pdf"))
    Plots.savefig(Plots.plot(RoGs; leg=false),
        joinpath(figure_path, "RoG.pdf"))

    for (snap, i) in zip(eachcol(snapshots), snapshot_inds)
        Plots.savefig(plot_chain(snap), joinpath(figure_path, "snapshot$i.pdf"))
    end
end

function run_2_1_7a()
    figure_path = "figures/2.1.7/a"

    interaction_matrix = make_interaction_energy_matrix()

    n = 100
    monomer_types = rand(1:20, n)

    chain = make_linear_2d_chain(n)
    coord_map = make_coord_map(chain)

    temperatures = [10, 8, 6, 4, 2, 1]
    sweeps = 100n
    interfaces = [i * sweeps for i in 1:length(temperatures)-1]

    energies, e2e_dists, RoGs = simulate(
        chain, coord_map, interaction_matrix, monomer_types, 10, 0
    )

    for temperature in temperatures
        energies, e2e_dists, RoGs = @time simulate(
            chain, coord_map, interaction_matrix,
            monomer_types, temperature, sweeps;
            energies=energies, e2e_dists=e2e_dists, RoGs=RoGs
        )
    end

    max_points = 1000
    spacing = cld(length(energies), max_points)
    x_axis = 1:spacing:length(energies)

    Plots.savefig(plot_temperature_interfaces!(
            Plots.plot(x_axis, (@view energies[x_axis]); leg=false),
            interfaces),
        joinpath(figure_path, "energy$n.pdf"))
    Plots.savefig(plot_temperature_interfaces!(
            Plots.plot(x_axis, (@view e2e_dists[x_axis]); leg=false),
            interfaces),
        joinpath(figure_path, "e2e$n.pdf"))
    Plots.savefig(plot_temperature_interfaces!(
            Plots.plot(x_axis, (@view RoGs[x_axis]); leg=false),
            interfaces),
        joinpath(figure_path, "RoG$n.pdf"))
end

function run_2_1_7b()
    function make_data(interaction_matrix, monomer_types, chain, coord_map,
        temperatures, sweeps_eq, sweeps_mean, sweeps_init)
        energy_t = Float64[]
        e2e_t = Float64[]
        RoG_t = Float64[]

        simulate_mean(
            chain, coord_map, interaction_matrix,
            monomer_types, temperatures[1], sweeps_init, 0
        )

        for temperature in temperatures
            energy_mean, e2e_mean, RoG_mean = simulate_mean(
                chain, coord_map, interaction_matrix,
                monomer_types, temperature, sweeps_eq, sweeps_mean
            )

            push!(energy_t, energy_mean)
            push!(e2e_t, e2e_mean)
            push!(RoG_t, RoG_mean)
        end

        energy_t, e2e_t, RoG_t
    end

    figure_path = "figures/2.1.7/b"

    interaction_matrix = make_interaction_energy_matrix()

    # n = 15  => init = 10_000n,  eq = 100n,  mean = 10_000n
    # n = 50  => init = 100_000n, eq = 1000n, mean = 100_000n
    # n = 100 => init = 100_000n, eq = 1000n, mean = 100_000n

    n = 100
    monomer_types = rand(1:20, n)

    # nth = Threads.nthreads()
    nth = 96

    chains = [make_linear_2d_chain(n) for _ in 1:nth]
    coord_maps = [make_coord_map(chain) for chain in chains]

    temperatures = range(20, 0.1, 100)
    sweeps_eq = 1000n
    sweeps_mean = 100_000n
    sweeps_init = 100_000n

    data = @time pmap(zip(chains, coord_maps)) do (chain, coord_map)
        @inbounds make_data(interaction_matrix, monomer_types,
            chain, coord_map, temperatures, sweeps_eq, sweeps_mean, sweeps_init)
    end

    energy_thread = [x for (x, _, _) in data]
    e2e_thread = [x for (_, x, _) in data]
    RoG_thread = [x for (_, _, x) in data]

    energy_t = mean(energy_thread)
    e2e_t = mean(e2e_thread)
    RoG_t = mean(RoG_thread)

    Plots.savefig(Plots.plot(temperatures, energy_t;
            leg=false, xlabel="T", ylabel="Energy"),
        joinpath(figure_path, "energy$n.pdf"))
    Plots.savefig(Plots.plot(temperatures, e2e_t;
            leg=false, xlabel="T", ylabel="E2E distance"),
        joinpath(figure_path, "e2e$n.pdf"))
    Plots.savefig(Plots.plot(temperatures, RoG_t;
            leg=false, xlabel="T", ylabel="RoG"),
        joinpath(figure_path, "RoG$n.pdf"))
end

function run_2_1_8()
    figure_path = "figures/2.1.8"

    interaction_matrix = make_interaction_energy_matrix()

    n = 30
    monomer_types = rand(1:20, n)

    sweeps = 1000n
    log_interval = sweeps ÷ 1000
    T_start = 10
    T_end = 0.5

    x_axis = 0:log_interval:sweeps

    for i in 1:2
        chain = make_linear_2d_chain(n)
        coord_map = make_coord_map(chain)

        energies = @time simulate_annealing(
            chain, coord_map, interaction_matrix, monomer_types,
            sweeps, 1, 1, log_interval
        )

        Plots.savefig(Plots.plot(x_axis, energies; leg=false),
            joinpath(figure_path, "energy_$i.pdf"))

        Plots.savefig(plot_chain(chain),
            joinpath(figure_path, "chain_$i.pdf"))
    end

    chain = make_linear_2d_chain(n)
    coord_map = make_coord_map(chain)

    energies = @time simulate_annealing(
        chain, coord_map, interaction_matrix, monomer_types,
        sweeps, T_start, T_end, log_interval
    )

    Plots.savefig(Plots.plot(x_axis, energies; leg=false),
        joinpath(figure_path, "energy_SA.pdf"))

    Plots.savefig(plot_chain(chain),
        joinpath(figure_path, "chain_SA.pdf"))
end

function run_2_1_9()
    figure_path = "figures/2.1.9"

    for p in (0.1, 0.5, 0.9)
        interaction_matrix = make_repulsive_interaction_energy_matrix(p)

        Plots.savefig(Plots.heatmap(interaction_matrix[:, end:-1:1]),
            joinpath(figure_path, "matrix_$p.pdf"))

        n = 50
        monomer_types = rand(1:20, n)

        sweeps = 10000n
        T_start = 10
        T_end = 0.5

        chain = make_linear_2d_chain(n)
        coord_map = make_coord_map(chain)

        @time simulate_annealing(
            chain, coord_map, interaction_matrix, monomer_types,
            sweeps, T_start, T_end
        )

        Plots.savefig(plot_chain(chain),
            joinpath(figure_path, "chain_$p.pdf"))
    end
end
