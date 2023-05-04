# 3d versions of utility functions to rotate a vector 90 degrees around an axis
rot3d(r, (x, y, z),) = (
    r .+ (x, -z, y),
    r .+ (x, z, -y),
    r .+ (-z, y, x),
    r .+ (z, y, -x),
    r .+ (-y, x, z),
    r .+ (y, -x, z),
)

# This function will check if moving monomer number i is possible, and return
# the destination coordinate if it is.
function is_transition_possible(chain::Vector{NTuple{3,Int}},
    coord_map, i)
    if i == 1 || i == length(chain)
        m1, m2 = if i == 1
            chain[1], chain[2]
        else
            chain[end], chain[end-1]
        end

        dests = @inline rot3d(m2, m1 .- m2)

        possible = ((!haskey(coord_map, dest) for dest in dests)...,)

        amt_possible = count(possible)

        if amt_possible > 0
            choice = rand(1:amt_possible)
            choice_dest = (0, 0, 0)

            i = 0
            for (dest, b) in zip(dests, possible)
                i += b
                if i == choice
                    choice_dest = dest
                    break
                end
            end

            return (true, choice_dest)
        end
    else
        # Short hand to not have to write out 'chain[...]' too many times
        m0, m1, m2 = chain[i-1], chain[i], chain[i+1]

        # We need to be on a corner, which means that the vector from the
        # monomer before to the monomer after needs to be different in
        # two coordinates.
        if count(m2 .!= m0) == 2
            # Observing that the new coordinate of the monomer is obtained from
            # moving it towards both its bonded monomers.
            dest_coord = m0 .+ m2 .- m1
            return (!haskey(coord_map, dest_coord), dest_coord)
        end
    end

    (false, (0, 0))
end
