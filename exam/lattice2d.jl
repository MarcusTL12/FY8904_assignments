# This is for initializing a 2d linear "protein" chain starting at the origin
# and reaching towards positive x
# Calling this for n = 4 will return the list
# [(0, 0), (1, 0), (2, 0), (3, 0)]
function make_linear_2d_chain(n)
    [(i, 0) for i in 0:(n-1)]
end

# This is for initializing a 2d perfectly folded protein starting at the origin
# and moving upwards 'w' units, then right one, then down another 'w' units.
# n = 20, and w = 5 gives the following chain:
# ╔╗╔╗║
# ║║║║║
# ║║║║║
# ║║║║║
# ║╚╝╚╝
function make_folded_2d_chain(n, w)
    coord = (0, 0)

    chain = NTuple{2,Int}[]

    dir = 1

    while length(chain) < n
        for i in 1:w
            push!(chain, coord)
            if i < w
                coord = coord .+ (0, dir)
            end

            if length(chain) >= n
                break
            end
        end

        coord = coord .+ (1, 0)

        dir = -dir
    end

    chain
end

# This translates the chain into a dictionary (Hashmap) from integer coordinates
# to the index in the chain. This makes for a both space and time efficient
# lookup table for whether a given coordinate is occupied, and by which one.
#
# This function is written very generally and will work for making a Dict
# from elements to indices for any array, most notably also for 3d coordinates
function make_coord_map(chain)
    Dict(coord => i for (i, coord) in enumerate(chain))
end

# This takes in the previous neighbour_list (or an empty one), and updates
# the list to the given chain. This is to omit the allocation of a new list.
function make_neighbour_list!(neighbour_list, chain, coord_map)
    # Remove all elements from list
    empty!(neighbour_list)

    # Loop over each link in the chain, and then each of the four nearest
    # neighbours.
    for (i, coord) in enumerate(chain)
        for dir in ((1, 0), (-1, 0), (0, 1), (0, -1))
            neighbour_coord = coord .+ dir

            # get the index in the chain (or 0 if it does not exist)
            # and then only add a neighbour pair if the neighbour index
            # is more than one above the current monomer. This means that it
            # is later in the chain (So that we only store one element per
            # neighbour pair), and they are not covalently bonded as they
            # are not adjacent in the chain.
            neighbour_ind = get(coord_map, neighbour_coord, 0)
            if neighbour_ind > i + 1
                push!(neighbour_list, (i, neighbour_ind))
            end
        end
    end

    neighbour_list
end

# Same as function above, but makes a new list
function make_neighbour_list(chain, coord_map)
    make_neighbour_list!(NTuple{2,Int}[], chain, coord_map)
end

# Calculate the total non-covalent energy from the neighbour_list
function calculate_energy(interaction_matrix, monomer_types, neighbour_list)
    sum(interaction_matrix[monomer_types[i], monomer_types[j]]
        for (i, j) in neighbour_list)
end

# Calculate the total energy directly without having to save the neighbour_list.
# The code here is almost identical to the make_neighbour_list! function
# now just adding to the total energy instead of adding elements to the
# neighbour_list.
#
# The explicit type for the chain parameter is to avoid function name
# collisions with the 3d version
function calculate_energy_direct(interaction_matrix, monomer_types,
    chain::AbstractArray{NTuple{2,Int}}, coord_map)

    energy = 0.0

    # Loop over each link in the chain, and then each of the four nearest
    # neighbours.
    for (i, coord) in enumerate(chain)
        for dir in ((1, 0), (-1, 0), (0, 1), (0, -1))
            neighbour_coord = coord .+ dir

            # get the index in the chain (or 0 if it does not exist)
            # and then only add a neighbour pair if the neighbour index
            # is more than one above the current monomer. This means that it
            # is later in the chain (So that we only store one element per
            # neighbour pair), and they are not covalently bonded as they
            # are not adjacent in the chain.
            neighbour_ind = get(coord_map, neighbour_coord, 0)
            if neighbour_ind > i + 1
                energy += interaction_matrix[monomer_types[i],
                    monomer_types[neighbour_ind]]
            end
        end
    end

    energy
end

# Calculate the energy contribution from monomer number i only. This is very
# useful for calculating the ΔE values during MC moves.
function get_energy_contrib_i(interaction_matrix, monomer_types,
    chain::AbstractArray{NTuple{2,Int}}, coord_map, i)

    coord = chain[i]

    energy = 0.0

    for dir in ((1, 0), (-1, 0), (0, 1), (0, -1))
        neighbour_coord = coord .+ dir

        j = get(coord_map, neighbour_coord, 0)
        if j > 0 && abs(i - j) > 1
            energy += interaction_matrix[monomer_types[i], monomer_types[j]]
        end
    end

    energy
end

# Function to calculate the euclidian end to end distance of a chain.
# Written generically enough to work for any number of dimensions.
function calculate_end2end_dist(chain)
    hypot((chain[1] .- chain[end])...)
end

# Calculates the center of "mass" of a chain assuming all amino acids weigh the
# same. This is exlusively used for calculating the radius of gyration.
function calculate_CoM(chain)
    # Initializing the (0.0, 0.0) tuple like this makes the function generic
    # enough to work for 3d as well
    CoM = chain[1] .* 0.0

    for coord in chain
        CoM = CoM .+ coord
    end

    CoM ./ length(chain)
end

# This calculates the radius of gyration of a chain
# The code is written generally enough to work for 3d,
# however in 3d it does not calculate the "classical" RoG, but the IUPAC version
# of the RoG for polymers. This is quite fitting for the purpose, but it is
# purely coincidence.
function calculate_RoG(chain)
    CoM = calculate_CoM(chain)

    r2 = 0.0

    for coord in chain, q in coord .- CoM
        r2 += q^2
    end

    √(r2 / length(chain))
end

# Returns a symmetric 20x20 matrix with elements uniformly distributed
# between -4 and -2
function make_interaction_energy_matrix()
    Symmetric([rand() * 2 - 4 for _ in 1:20, _ in 1:20])
end

# This works like the one above but takes in the probability for the elements to
# flip sign.
function make_repulsive_interaction_energy_matrix(prob)
    Symmetric([(rand() * 2 - 4) * ((rand() > prob) * 2 - 1)
               for _ in 1:20, _ in 1:20])
end

