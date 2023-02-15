using LinearAlgebra

using CairoMakie
using GLMakie

function make_arrow_figure(lattice_points, spin)
    f = Figure(resolution=(1000, 1000))

    ax = Axis3(f[1, 1], viewmode=:fit, aspect=:data, perspectiveness=0.5)

    arrows!(ax,
        (@view lattice_points[:, 1]),
        (@view lattice_points[:, 2]),
        (@view lattice_points[:, 3]),
        (@view spin[:, 1]),
        (@view spin[:, 2]),
        (@view spin[:, 3]);
        lengthscale=0.2, linewidth=0.05, arrowsize=0.1, color=:gray
    )

    f
end

function test_makie(n)
    lattice_points = zeros(n, n, n, 3)

    spin = randn(n, n, n, 3)

    for i in 1:n, j in 1:n, k in 1:n
        lattice_points[i, j, k, 1] = Float64(i)
        lattice_points[i, j, k, 2] = Float64(j)
        lattice_points[i, j, k, 3] = Float64(k)

        v = @view spin[i, j, k, :]
        v ./= norm(v)
    end

    lattice_points = reshape(lattice_points, n^3, 3)
    spin = reshape(spin, n^3, 3)

    make_arrow_figure(lattice_points, spin)
end

function test_anim()
    f = Figure(resolution=(1000, 1000))

    ax = Axis(f[1, 1])

    n = 1000
    nt = 1000

    xs = range(0, 2π, n)

    data = zeros(n, nt)

    for i in 1:nt-1
        data[:, i+1] .= (@view data[:, i]) .+ randn(n) .* 0.01
    end

    ys = Observable(@view data[:, 1])

    lines!(ax, xs, ys)

    record(f, "tmp.mp4", 1:nt, framerate=60) do i
        ys[] = @view data[:, i]
    end
end
