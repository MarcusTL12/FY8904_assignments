
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

function make_x_neighbour_ranges(lattice)
    x_neighbour_ranges = UnitRange{Int64}[]

    start = -1

    for i in 1:length(lattice)-1
        if start == -1
            if lattice[i] > 0
                start = lattice[i]
            end
        elseif lattice[i] <= 0
            stop = lattice[i - 1]
            push!(x_neighbour_ranges, start:stop)
            start = -1
        end
    end

    x_neighbour_ranges
end
