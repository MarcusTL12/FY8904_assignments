
function test_3d()
    interaction_matrix = make_interaction_energy_matrix()

    n = 30
    monomer_types = rand(1:20, n)

    chain = make_linear_3d_chain(n)
    coord_map = make_coord_map(chain)

    @time simulate_annealing(chain, coord_map, interaction_matrix,
        monomer_types, 10000n, 10, 0.5)

    plot_chain(chain)
end

function run_2_2_2()
    figure_path = "figures/2.2.2"

    interaction_matrix = make_interaction_energy_matrix()

    n = 15
    monomer_types = rand(1:20, n)

    chain = make_linear_3d_chain(n)
    coord_map = make_coord_map(chain)

    snapshot_inds = (1, 10, 100)

    _, _, _, snapshots = @time simulate(
        chain, coord_map, interaction_matrix, monomer_types, 10, 100,
        snapshot_inds
    )

    Plots.plotly()
    for (snap, i) in zip(eachcol(snapshots), snapshot_inds)
        Plots.savefig(plot_chain(snap),
            joinpath(figure_path, "snapshot$i.html"))
    end
    Plots.gr()
end

function run_2_2_3()
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

    figure_path = "figures/2.2.3/"

    interaction_matrix = make_interaction_energy_matrix()

    # n = 15  => init = 10_000n, eq = 100n, mean = 10_000n
    # n = 50  => inti = 10_000n, eq = 100n, mean = 10_000n
    # n = 100 => init = 10_000n, eq = 100n, mean = 10_000n

    n = 15
    monomer_types = rand(1:20, n)

    # nth = Threads.nthreads()
    nth = 96

    chains = [make_linear_3d_chain(n) for _ in 1:nth]
    coord_maps = [make_coord_map(chain) for chain in chains]

    temperatures = range(20, 0.1, 100)
    sweeps_eq = 100n
    sweeps_mean = 10_000n
    sweeps_init = 10_000n

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

function run_annealing_ljp()
    interaction_matrix = make_interaction_energy_matrix()

    n = 50
    monomer_types = rand(1:20, n)

    chain = make_linear_3d_chain(n)
    coord_map = make_coord_map(chain)

    @show calculate_ljp_energy(interaction_matrix, monomer_types, chain)

    @time simulate_ljp_annealing(chain, coord_map, interaction_matrix,
        monomer_types, 10000n, 10, 0.5)

    plot_chain(chain)
end

function run_ljp_phase_diagram()
    function make_data(interaction_matrix, monomer_types, chain, coord_map,
        temperatures, sweeps_eq, sweeps_mean, sweeps_init)
        energy_t = Float64[]
        e2e_t = Float64[]
        RoG_t = Float64[]

        simulate_ljp_mean(
            chain, coord_map, interaction_matrix,
            monomer_types, temperatures[1], sweeps_init, 0
        )

        for temperature in temperatures
            energy_mean, e2e_mean, RoG_mean = simulate_ljp_mean(
                chain, coord_map, interaction_matrix,
                monomer_types, temperature, sweeps_eq, sweeps_mean
            )

            push!(energy_t, energy_mean)
            push!(e2e_t, e2e_mean)
            push!(RoG_t, RoG_mean)
        end

        energy_t, e2e_t, RoG_t
    end

    figure_path = "figures/ljp/"

    interaction_matrix = make_interaction_energy_matrix()

    # n = 15  => init = 10_000n, eq = 1000n, mean = 10_000n
    # n = 50  => inti = 1000n, eq = 100n, mean = 1000n
    # n = 100 => init = 1000n, eq = 100n, mean = 1000n

    n = 100
    monomer_types = rand(1:20, n)

    # nth = Threads.nthreads()
    nth = 96

    chains = [make_linear_3d_chain(n) for _ in 1:nth]
    coord_maps = [make_coord_map(chain) for chain in chains]

    temperatures = range(30, 0.1, 100)
    sweeps_eq = 100n
    sweeps_mean = 1000n
    sweeps_init = 1000n

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
