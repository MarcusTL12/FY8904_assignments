
# Simple utility functions to rotate a vector 90 degrees left or right
rotl((x, y),) = (-y, x)
rotr((x, y),) = (y, -x)

# This function will check if moving monomer number i is possible, and return
# the destination coordinate if it is.
function is_transition_possible(chain::AbstractArray{NTuple{2,Int}},
    coord_map, i)
    if i == 1 || i == length(chain)
        dest1, dest2 = if i == 1
            chain[2] .+ rotl(chain[1] .- chain[2]),
            chain[2] .+ rotr(chain[1] .- chain[2])
        else
            chain[end-1] .+ rotl(chain[end] .- chain[end-1]),
            chain[end-1] .+ rotr(chain[end] .- chain[end-1])
        end

        possible1 = !haskey(coord_map, dest1)
        possible2 = !haskey(coord_map, dest2)

        if possible1 && possible2
            # If both are possible, pick one at random
            return (true, rand((dest1, dest2),))
        elseif possible1
            return (true, dest1)
        elseif possible2
            return (true, dest2)
        end
    else
        # Short hand to not have to write out 'chain[...]' too many times
        m0, m1, m2 = chain[i-1], chain[i], chain[i+1]

        # Semi black magic to check whether monomer i is on a corner
        # By forming the vector from monomer i - 1 and i + 1, this needs to
        # have a non-zero x and y component for i to be on a corner.
        # We check this by multiplying the x and y component and checking
        # if the result is 0
        if prod(m2 .- m0) != 0
            # Observing that the new coordinate of the monomer is obtained from
            # moving it towards both its bonded monomers.
            dest_coord = m0 .+ m2 .- m1
            return (!haskey(coord_map, dest_coord), dest_coord)
        end
    end

    (false, (0, 0))
end

function do_mc_draw!(chain, coord_map, interaction_matrix, monomer_types,
    prev_energy, temperature)

    # Draw one of the possible transitions at random.
    i = rand(1:length(chain))

    possible, dest_coord = is_transition_possible(chain, coord_map, i)

    if possible
        # Save this to revert if rejected
        prev_coord = chain[i]

        E1 = get_energy_contrib_i(interaction_matrix, monomer_types, chain,
            coord_map, i)

        # Update the structure
        chain[i] = dest_coord
        delete!(coord_map, prev_coord)
        coord_map[dest_coord] = i

        E2 = get_energy_contrib_i(interaction_matrix, monomer_types, chain,
            coord_map, i)

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
function do_mc_sweep!(chain, coord_map, interaction_matrix, monomer_types,
    energy, temperature)
    for _ in 1:length(chain)
        energy = do_mc_draw!(chain, coord_map, interaction_matrix,
            monomer_types, energy, temperature)
    end
    energy
end

function simulate(chain, coord_map, interaction_matrix, monomer_types,
    temperature, n, snapshot_inds=();
    energies=[calculate_energy_direct(interaction_matrix, monomer_types,
        chain, coord_map)],
    e2e_dists=[calculate_end2end_dist(chain)],
    RoGs=[calculate_RoG(chain)],
    snapshots=eltype(chain)[]
)

    energy = energies[end]

    for i in 1:n
        energy = do_mc_sweep!(chain, coord_map, interaction_matrix,
            monomer_types, energy, temperature)

        push!(energies, energy)
        push!(e2e_dists, calculate_end2end_dist(chain))
        push!(RoGs, calculate_RoG(chain))

        if i ∈ snapshot_inds
            append!(snapshots, chain)
        end
    end

    snapshots = reshape(snapshots, length(chain), length(snapshot_inds))

    energies, e2e_dists, RoGs, snapshots
end

function simulate_mean(chain, coord_map, interaction_matrix,
    monomer_types, temperature, n_eq, n_mean)

    energy = calculate_energy_direct(interaction_matrix, monomer_types,
        chain, coord_map)

    for _ in 1:n_eq
        energy = do_mc_sweep!(chain, coord_map, interaction_matrix,
            monomer_types, energy, temperature)
    end

    energy_mean = 0.0
    e2e_mean = 0.0
    RoG_mean = 0.0

    for _ in 1:n_mean
        energy = do_mc_sweep!(chain, coord_map, interaction_matrix,
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

function simulate_annealing(chain, coord_map, interaction_matrix,
    monomer_types, n, T_start, T_end, log_interval=0)

    energy = calculate_energy_direct(interaction_matrix, monomer_types,
        chain, coord_map)

    energies = [energy]

    for (i, temperature) in enumerate(range(T_start, T_end, n))
        energy = do_mc_sweep!(chain, coord_map, interaction_matrix,
            monomer_types, energy, temperature)

        if log_interval != 0 && i % log_interval == 0
            push!(energies, energy)
        end
    end

    energies
end
