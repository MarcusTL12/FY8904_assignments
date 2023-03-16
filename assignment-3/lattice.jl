
const outside = 0
const boundary = 1
const interior = 2

function setup_lattice(corners, n)
    minx, maxx = extrema(real, corners)
    miny, maxy = extrema(imag, corners)

    minx, maxx, miny, maxy

    xl = length(minx:maxx)
    yl = length(miny:maxy)

    xn = (xl - 1) * n + 1
    yn = (yl - 1) * n + 1

    lattice = zeros(Int32, xn, yn)

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

    while !isempty(queue)
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

function test_lattice(l, n)
    corners = make_koch_square(l)
    setup_lattice(corners, n)
end
