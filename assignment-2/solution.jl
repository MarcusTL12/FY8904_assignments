include("hamiltonian.jl")
include("simulation.jl")
include("visualization.jl")
using FFTW

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

function test_1d_chain(n)
    S = zeros(n, 1, 1, 3)
    # S = randn(n, 1, 1, 3)

    # randn!(@view S[:, :, 1, 1:2])
    # S .*= 0.2

    S[:, :, :, 3] .= 1.0
    # S[n÷2, 1, 1, 2] = 1.0
    # S[n÷2, 1, 1, 3] = 0.0

    normalize_spin!(S)

    lattice_points = zeros(n, 3)

    for x in 1:n
        lattice_points[x, 1] = Float64(x)
        lattice_points[x, 2] = 0.0
        lattice_points[x, 3] = 0.0
    end

    state = init_state(S)
    params = setup_params(
        10.0, 3.0, 1.0, (@SVector [0.0, 0.0, 0.0]), 1.0, 0.1
    )

    S_hist = @time simulate!(state, params, 1000, 1)

    @time visualize_spin_history_interactive(lattice_points, S_hist)
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

    # randn!(@view S[:, :, 1, 1:2])
    # S .*= 0.2

    # @. S[:, :, 1, 3] =
    #     √(1.0 - (@view S[:, :, 1, 2])^2 - (@view S[:, :, 1, 3])^2)

    # S[:, :, :, 3] .= 1.0
    # S[1, 1, 1, 1] = 1.0
    # S[1, 1, 1, 3] = 0.0

    normalize_spin!(S)

    lattice_points = zeros(nx, ny, nz, 3)

    for x in 1:nx, y in 1:ny, z in 1:nz
        lattice_points[x, y, z, 1] = Float64(x)
        lattice_points[x, y, z, 2] = Float64(y)
        lattice_points[x, y, z, 3] = Float64(z)
    end

    lattice_points = reshape(lattice_points, n, 3)

    state = init_state(S)
    params = setup_params(
        10.0, 6.0, 0.1, (@SVector [0.0, 0.0, 0.0]), 1.0, 0.1
    )

    S_hist = @time simulate!(state, params, 4000, 1)

    @time visualize_spin_history_interactive(lattice_points, S_hist)
end

function test_3d_box_nstaged(nx, ny, nz)
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
        -10.0, 0.0, 0.1, (@SVector [0.0, 0.0, 0.0]), 0.1, 0.1
    )

    S_hist = @time simulate!(state, params, 500, 10)

    for B in 10.0:10.0:200.0
        params = setup_params(
            -10.0, 0.0, 0.1, (@SVector [0.0, 0.0, B]), 0.1, 0.1
        )

        S_hist = @time simulate!(state, params, 100, 10, S_hist, false)
    end

    params = setup_params(
        -10.0, 0.0, 0.1, (@SVector [0.0, 0.0, 0.0]), 0.1, 0.1
    )

    S_hist = @time simulate!(state, params, 1000, 10, S_hist, false)

    @time visualize_spin_history_interactive(lattice_points, S_hist)
end

function test_1d_dispersion(n)
    S = zeros(n, 1, 1, 3)

    S[:, :, :, 3] .= 1.0

    # normalize_spin!(S)

    # lattice_points = zeros(n, 3)

    # for x in 1:n
    #     lattice_points[x, 1] = Float64(x)
    #     lattice_points[x, 2] = 0.0
    #     lattice_points[x, 3] = 0.0
    # end

    state = init_state(S)
    J = 10.0
    dz = 3.0
    params = setup_params(
        J, dz, 0.01, (@SVector [0.0, 0.0, 0.0]), 1.0, 0.1
    )

    n_steps = 30_000
    S_hist = @time simulate!(state, params, n_steps, 1)

    Sx = @view S_hist[:, 1, :]

    # h/2 = 2067.811956368851 meV * fs
    f_analytic(k) = (2dz + 2J * (1 - cos(2π * k))) / 2067.811956368851

    Sx_fft = @time fftshift(fft(Sx))
    k_fft = fftshift(fftfreq(n))
    f_fft = fftshift(fftfreq(n_steps + 1), 1)

    f, ax, _ = heatmap(k_fft, f_fft, norm.(Sx_fft))
    plot!(ax, k_fft, f_analytic)

    ax.xlabel[] = "k"
    ax.ylabel[] = "f"

    f
end

function test_1d_dispersion_antiferromagnet(n)
    S = zeros(n, 1, 1, 3)

    S[1:2:end, :, :, 3] .= 1.0
    S[2:2:end, :, :, 3] .= -1.0

    # normalize_spin!(S)

    # lattice_points = zeros(n, 3)

    # for x in 1:n
    #     lattice_points[x, 1] = Float64(x)
    #     lattice_points[x, 2] = 0.0
    #     lattice_points[x, 3] = 0.0
    # end

    state = init_state(S)
    J = -30.0
    dz = 6.0
    params = setup_params(
        J, dz, 0.01, (@SVector [0.0, 0.0, 0.05J]), 1.0, 0.1
    )

    n_steps = 30_000
    S_hist = @time simulate!(state, params, n_steps, 1)

    Sx = @view S_hist[:, 1, :]

    # h/2 = 2067.811956368851 meV * fs
    f_analytic(k) = (2dz + 2J * (1 - cos(2π * k))) / 2067.811956368851

    Sx_fft = @time fftshift(fft(Sx))
    k_fft = fftshift(fftfreq(n))
    f_fft = fftshift(fftfreq(n_steps + 1), 1)

    f, ax, _ = heatmap(k_fft, f_fft, norm.(Sx_fft))
    # plot!(ax, k_fft, f_analytic)

    ax.xlabel[] = "k"
    ax.ylabel[] = "f"

    f
end
