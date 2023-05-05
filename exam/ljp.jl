
# Defining the sigma value of the LJP such that the equilibrium distance
# is equal to 1
const σ = 2^(-1 / 6)

# Function to calculate the lennard jones potential efficently
function lennard_jones(r1, r2, ε)
    r = @inline hypot((r1 .- r2)...)

    rs = σ / r
    rs3 = rs * rs * rs
    rs6 = rs3 * rs3
    rs12 = rs6 * rs6

    4ε * (rs12 - rs6)
end

function calculate_ljp_energy(interaction_matrix, monomer_types, chain)
    energy = 0.0

    @inbounds for (i, coord_i) in enumerate(chain), j in i+1:length(chain)
        energy += lennard_jones(coord_i, chain[j],
            interaction_matrix[monomer_types[i], monomer_types[j]])
    end

    energy
end

function get_ljp_energy_contrib_i(interaction_matrix, monomer_types, chain, i)
    energy = 0.0

    coord_i = @inbounds chain[i]

    @inbounds for (j, coord_j) in enumerate(chain)
        energy += lennard_jones(coord_i, coord_j,
            interaction_matrix[monomer_types[i], monomer_types[j]])
    end

    energy
end

function do_mc_ljp_draw!(chain, coord_map, interaction_matrix, monomer_types,
    prev_energy, temperature)

    # Draw one of the possible transitions at random.
    i = rand(1:length(chain))

    possible, dest_coord = is_transition_possible(chain, coord_map, i)

    if possible
        # Save this to revert if rejected
        prev_coord = chain[i]

        E1 = get_ljp_energy_contrib_i(interaction_matrix, monomer_types,
            chain, i)

        # Update the structure
        chain[i] = dest_coord
        delete!(coord_map, prev_coord)
        coord_map[dest_coord] = i

        E2 = get_ljp_energy_contrib_i(interaction_matrix, monomer_types,
            chain, i)

        # Calculate the ΔE directly from the energy contribution of the
        # monomer that is being moved.
        ΔE = E2 - E1
        new_energy = prev_energy + ΔE

        # Metropolis acceptance probability
        p_accept = min(1.0, exp(-ΔE / temperature))

        # If rejected, we need to revert the changes to the polymer structure
        if rand() > p_accept
            chain[i] = prev_coord
            delete!(coord_map, dest_coord)
            coord_map[prev_coord] = i
            prev_energy
        else
            new_energy
        end
    else
        prev_energy
    end
end

# Just do N draws
function do_mc_ljp_sweep!(chain, coord_map, interaction_matrix, monomer_types,
    energy, temperature)
    for _ in 1:length(chain)
        energy = do_mc_ljp_draw!(chain, coord_map, interaction_matrix,
            monomer_types, energy, temperature)
    end
    energy
end

function simulate_ljp_mean(chain, coord_map, interaction_matrix,
    monomer_types, temperature, n_eq, n_mean)

    energy = calculate_ljp_energy(interaction_matrix, monomer_types, chain)

    for _ in 1:n_eq
        energy = do_mc_ljp_sweep!(chain, coord_map, interaction_matrix,
            monomer_types, energy, temperature)
    end

    energy_mean = 0.0
    e2e_mean = 0.0
    RoG_mean = 0.0

    for _ in 1:n_mean
        energy = do_mc_ljp_sweep!(chain, coord_map, interaction_matrix,
            monomer_types, energy, temperature)

        energy_mean += energy
        e2e_mean += calculate_end2end_dist(chain)
        RoG_mean += calculate_RoG(chain)
    end

    energy_mean /= n_mean
    e2e_mean /= n_mean
    RoG_mean /= n_mean

    energy_mean, e2e_mean, RoG_mean
end

function simulate_ljp_annealing(chain, coord_map, interaction_matrix,
    monomer_types, n, T_start, T_end, log_interval=0)

    energy = calculate_ljp_energy(interaction_matrix, monomer_types, chain)

    energies = [energy]

    for (i, temperature) in enumerate(range(T_start, T_end, n))
        energy = do_mc_ljp_sweep!(chain, coord_map, interaction_matrix,
            monomer_types, energy, temperature)

        if log_interval != 0 && i % log_interval == 0
            push!(energies, energy)
        end
    end

    energies
end
