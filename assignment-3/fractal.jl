using Plots

const koch_generator = (
    1 + 0im,
    0 + 1im,
    1 + 0im,
    0 - 1im,
    0 - 1im,
    1 + 0im,
    0 + 1im,
    1 + 0im
)

function generate_koch!(corners, dir, l)
    if l == 0
        push!(corners, last(corners) + dir)
    else
        for d in koch_generator
            generate_koch!(corners, d * dir, l - 1)
        end
        corners
    end
end

function make_koch_square(l)
    corners = [0 + 0im]
    generate_koch!(corners, 1 + 0im, l)
    generate_koch!(corners, 0 + 1im, l)
    generate_koch!(corners, -1 + 0im, l)
    generate_koch!(corners, 0 - 1im, l)
    corners
end

function plot_curve(corners)
    n_corners = length(corners)
    corners = reshape(reinterpret(Int, corners), 2, n_corners)

    xs = @view corners[1, :]
    ys = @view corners[2, :]

    plot(xs, ys, leg=false, ratio=1)
end

function plot_curve_scaled(corners)
    n_corners = length(corners)
    corners = reshape(reinterpret(Int, corners), 2, n_corners)

    xs = @view corners[1, :]
    ys = @view corners[2, :]

    L = maximum(xs)
    xs /= L
    ys /= L

    plot!(xs, ys)
end

function test_koch(l)
    plot(leg=false, ratio=1)

    for ml in 0:l
        corners = [0 + 0im]
        generate_koch!(corners, 1 + 0im, ml)
        plot_curve_scaled(corners)
    end

    plot!()
end
