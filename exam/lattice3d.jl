# This is for initializing a 3d linear "protein" chain starting at the origin
# and reaching towards positive x
# Calling this for n = 4 will return the list
# [(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)]
function make_linear_3d_chain(n)
    [(i, 0, 0) for i in 0:(n-1)]
end

# The make_coord_map function will work nicely with 3d coordinates as is
# so there is no need rewriting a 3d version

# I ended up not using the neighbour list implementation of the 2d chains,
# so I will not be implementing a 3d version of those.

# 3d version of the function to calculate the polymer energy
# Only difference being that here we loop over all 6 nearest neighbour coords
function calculate_energy_direct(interaction_matrix, monomer_types,
    chain::AbstractArray{NTuple{3,Int}}, coord_map)

    energy = 0.0

    # Loop over each link in the chain, and then each of the four nearest
    # neighbours.
    for (i, coord) in enumerate(chain)
        for dir in (
            (1, 0, 0), (-1, 0, 0),
            (0, 1, 0), (0, -1, 0),
            (0, 0, 1), (0, 0, -1)
        )
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

# 3d version with again the only difference being looping over
# 6 nearest neighbour coordinates
function get_energy_contrib_i(interaction_matrix, monomer_types,
    chain::AbstractArray{NTuple{3,Int}}, coord_map, i)

    coord = chain[i]

    energy = 0.0

    for dir in (
        (1, 0, 0), (-1, 0, 0),
        (0, 1, 0), (0, -1, 0),
        (0, 0, 1), (0, 0, -1)
    )
        neighbour_coord = coord .+ dir

        j = get(coord_map, neighbour_coord, 0)
        if j > 0 && abs(i - j) > 1
            energy += interaction_matrix[monomer_types[i], monomer_types[j]]
        end
    end

    energy
end
