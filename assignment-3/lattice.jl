
const outside = -2
const boundary = -1
const interior = 0

function setup_lattice(corners, n)
    minx, maxx = extrema(real, corners)
    miny, maxy = extrema(imag, corners)

    minx, maxx, miny, maxy

    xl = length(minx:maxx)
    yl = length(miny:maxy)

    xn = (xl - 1) * n + 1
    yn = (yl - 1) * n + 1

    lattice = fill(outside, xn, yn)

    x = (real(first(corners)) - minx) * n + 1
    y = (imag(first(corners)) - miny) * n + 1

    lattice[x, y] = 1

    for i in 1:length(corners)-1
        d = corners[i+1] - corners[i]
        for _ in 1:n
            x += real(d)
            y += imag(d)
            lattice[x, y] = boundary
        end
    end

    fill_interior!(lattice)
    index_lattice!(lattice)

    lattice
end

function fill_interior!(lattice)
    dirs = (
        (1, 0),
        (-1, 0),
        (0, 1),
        (0, -1),
    )

    queue = Deque{NTuple{2,Int}}()
    start = cld.(size(lattice), 2)
    lattice[start...] = interior
    push!(queue, start)

    @inbounds while !isempty(queue)
        p = popfirst!(queue)
        for d in dirs
            np = p .+ d
            if lattice[np...] == outside
                lattice[np...] = interior
                push!(queue, np)
            end
        end
    end
end

function index_lattice!(lattice)
    j = 1

    for (i, x) in enumerate(lattice)
        if x == interior
            lattice[i] = j
            j += 1
        end
    end
end

function make_lattice(l, n)
    corners = make_koch_square(l)
    setup_lattice(corners, n)
end

function get_h(l, n)
    1 / (4^l * n)
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
