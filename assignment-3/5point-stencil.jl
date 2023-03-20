function make_full_matrix(lattice, h)
    n = length(lattice)

    A = spzeros(n, n)

    ij_to_lin(i, j) = (j - 1) * size(lattice, 1) + i

    hinv2 = 1 / h^2

    for i in axes(A, 1)
        if lattice[i] > 0
            A[i, i] = 4hinv2
        end
    end

    for j in 1:size(lattice, 2)-1
        for i in 1:size(lattice, 1)-1
            if lattice[i, j] > 0 && lattice[i+1, j] > 0
                I = ij_to_lin(i, j)
                J = ij_to_lin(i + 1, j)
                A[I, J] = -hinv2
                A[J, I] = -hinv2
            end

            if lattice[i, j] > 0 && lattice[i, j+1] > 0
                I = ij_to_lin(i, j)
                J = ij_to_lin(i, j + 1)
                A[I, J] = -hinv2
                A[J, I] = -hinv2
            end
        end
    end

    A
end

function reduce_full_matrix(A, lattice)
    in_inds = [i for (i, x) in enumerate(lattice) if x > 0]

    A[in_inds, in_inds]
end

function make_reduced_matrix(lattice, h)
    hinv2 = 1 / h^2

    w = size(lattice, 1)

    d0 = [4hinv2 for _ in 1:length(lattice)]
    d1 = [-hinv2 for _ in 1:length(lattice)-1]
    dw = [-hinv2 for _ in 1:length(lattice)-w]

    A = spdiagm(0 => d0, 1 => d1, -1 => d1, w => dw, -w => dw)

    reduce_full_matrix(A, lattice)
end

# Returns an array of ranges of lattice points that are adjacent in the x
# direction
function make_x_neighbour_ranges(lattice)
    xnr = UnitRange{Int64}[]

    start = -1

    for (i, x) in enumerate(lattice)
        if start == -1
            if x > 0
                start = x
            end
        elseif x <= 0
            stop = lattice[i-1]
            push!(xnr, start:stop)
            start = -1
        end
    end

    xnr
end

# Returns an array of 3 ints of the structure (top_start, bottom_start, length)
function make_y_neighbour_ranges(lattice)
    ynr = NTuple{3,Int}[]

    top_start = -1
    bottom_start = -1

    for j in 2:size(lattice, 2)-1
        for i in axes(lattice, 1)
            x = lattice[i, j]
            y = lattice[i, j+1]
            if top_start == -1
                if x > 0 && y > 0
                    top_start = x
                    bottom_start = y
                end
            elseif x <= 0 || y <= 0
                top_stop = lattice[i-1, j]
                l = length(top_start:top_stop)
                push!(ynr, (top_start, bottom_start, l))
                top_start = -1
                bottom_start = -1
            end
        end
    end

    ynr
end


# Todo: add functions to make similar neighbour ranges for more distant
# than adjacent for higher order finite difference stencils.

function range_chunks(r, n)
    @inbounds begin
        rv = @inline r[1:n:end-(n-1)]
        rr = @inline r[length(rv)*n+1:end]
        rv, rr
    end
end

# Computes
# -∇u = A * u
# in place. A is represented by xnr, ynr and hinv
#
# n∇u = -∇u
# hinv2 = 1/h^2
#
# -∇u[i,j] = (4u[i,j] - u[i+1,j] - u[i-1,j] - u[i,j+1] - u[i,j-1]) / h^2
function compute_5p_laplacian!(n∇u, u, xnr, ynr, hinv2)
    # 4u[i, j] contribution
    @inbounds for i in eachindex(n∇u, u)
        n∇u[i] = u[i] * 4hinv2
    end

    # - u[i+1,j] - u[i-1,j] contribution
    @inbounds @fastmath for r in xnr
        @simd for i in r[1:end-1]
            n∇u[i] -= hinv2 * u[i+1]
        end
        @simd for i in r[1:end-1]
            n∇u[i+1] -= hinv2 * u[i]
        end
    end

    # - u[i,j+1] - u[i,j-1] contribution
    @inbounds @fastmath for (ts, bs, l) in ynr
        boff = bs - ts
        @simd for i in ts:ts+l-1
            j = i + boff
            n∇u[i] -= hinv2 * u[j]
            n∇u[j] -= hinv2 * u[i]
        end
    end
end

function make_linmap(lattice, l, n)
    h = get_h(l, n)

    xnr = make_x_neighbour_ranges(lattice)
    ynr = make_y_neighbour_ranges(lattice)
    hinv2 = 1 / h^2

    LinearMap(maximum(lattice);
        issymmetric=true, ismutating=true, isposdef=true) do n∇u, u
        compute_5p_laplacian!(n∇u, u, xnr, ynr, hinv2)
        n∇u
    end
end
