include("hamiltonian.jl")
include("simulation.jl")
include("visualization.jl")

function test_single_spin()
    Sx0 = randn() * 0.2
    Sy0 = randn() * 0.2
    Sz0 = √(1.0 - Sx0^2 - Sy0^2)
    S = zeros(1, 1, 1, 3)
    S[1, 1, 1, 1] = Sx0
    S[1, 1, 1, 2] = Sy0
    S[1, 1, 1, 3] = Sz0

    lattice_points = zeros(1, 3)

    state = init_state(S)
    params = setup_params(
        0.0, 1.0, 0.0, (@SVector [0.0, 0.0, 0.0]), 1.0, 0.0
    )

    S_hist = @time simulate!(state, params, 500, 10)

    @time animate_spin_history(lattice_points, S_hist)
end

function test_5x5_grid()
    S = zeros(5, 5, 1, 3)

    # randn!(@view S[:, :, 1, 1:2])
    # S .*= 0.2

    # @. S[:, :, 1, 3] =
    #     √(1.0 - (@view S[:, :, 1, 2])^2 - (@view S[:, :, 1, 3])^2)

    S[:, :, 1, 1] .= 1.0

    lattice_points = zeros(5, 5, 3)

    for x in 1:5, y in 1:5
        lattice_points[x, y, 1] = Float64(x)
        lattice_points[x, y, 2] = Float64(y)
    end

    lattice_points = reshape(lattice_points, 25, 3)

    state = init_state(S)
    params = setup_params(
        -10.0, 1.0, 1.0, (@SVector [0.0, 0.0, 0.0]), 1.0, 0.1
    )

    S_hist = @time simulate!(state, params, 1000, 1)

    @time animate_spin_history(lattice_points, S_hist)
end

function test_nxnxn_box(n)
    S = zeros(n, n, n, 3)

    # randn!(@view S[:, :, 1, 1:2])
    # S .*= 0.2

    # @. S[:, :, 1, 3] =
    #     √(1.0 - (@view S[:, :, 1, 2])^2 - (@view S[:, :, 1, 3])^2)

    S[:, :, :, 1] .= 1.0

    lattice_points = zeros(n, n, n, 3)

    for x in 1:n, y in 1:n, z in 1:n
        lattice_points[x, y, z, 1] = Float64(x)
        lattice_points[x, y, z, 2] = Float64(y)
        lattice_points[x, y, z, 3] = Float64(z)
    end

    lattice_points = reshape(lattice_points, n^3, 3)

    state = init_state(S)
    params = setup_params(
        -2.0, 1.0, 10.0, (@SVector [0.0, 0.0, 0.0]), 1.0, 0.1
    )

    S_hist = @time simulate!(state, params, 1000, 5)

    @time visualize_spin_history_interactive(lattice_points, S_hist)
end

function test_3d_box(nx, ny, nz)
    n = nx * ny * nz
    # S = zeros(nx, ny, nz, 3)
    S = randn(nx, ny, nz, 3)
    normalize_spin!(S)

    # randn!(@view S[:, :, 1, 1:2])
    # S .*= 0.2

    # @. S[:, :, 1, 3] =
    #     √(1.0 - (@view S[:, :, 1, 2])^2 - (@view S[:, :, 1, 3])^2)

    # S[:, :, :, 3] .= 1.0

    lattice_points = zeros(nx, ny, nz, 3)

    for x in 1:nx, y in 1:ny, z in 1:nz
        lattice_points[x, y, z, 1] = Float64(x)
        lattice_points[x, y, z, 2] = Float64(y)
        lattice_points[x, y, z, 3] = Float64(z)
    end

    lattice_points = reshape(lattice_points, n, 3)

    state = init_state(S)
    params = setup_params(
        10.0, 0.0, 1.0, (@SVector [0.0, 0.0, 0.0]), 0.1, 0.1
    )

    S_hist = @time simulate!(state, params, 500, 10)

    params = setup_params(
        10.0, 0.0, 1.0, (@SVector [0.0, 0.0, 100.0]), 0.1, 0.1
    )

    S_hist = @time simulate!(state, params, 500, 10, S_hist, false)

    @time visualize_spin_history_interactive(lattice_points, S_hist)
end
