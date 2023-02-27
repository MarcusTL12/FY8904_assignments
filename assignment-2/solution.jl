using LinearAlgebra
using ForwardDiff
using SIMD
using LoopVectorization
using Statistics
using StatsBase
using FFTW
FFTW.set_num_threads(Threads.nthreads())

using Random
using StaticArrays
# using GLMakie
# GLMakie.activate!(title="Assignment-2", framerate=60.0)
using CairoMakie

include("hamiltonian.jl")
include("simulation.jl")
include("visualization.jl")
include("magnetization.jl")

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

function test_1d_dispersion(n, gain=1.0, cutoff=0.0)
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
    Δt = 1.0
    params = setup_params(
        J, dz, 0.1, (@SVector [0.0, 0.0, 0.0]), Δt, 0.1
    )

    n_steps = ceil(Int, 30n / Δt)
    S_hist = @time simulate!(state, params, n_steps, 1)

    Sx = @view S_hist[:, 1, :]

    # h/2 = 2067.811956368851 meV * fs
    f_analytic(k) = (2dz + 2J * (1 - cos(2π * k))) / 2.067811956368851

    Sx_fft = @time fftshift(rfft(Sx')', 1)
    k_fft = fftshift(fftfreq(n))
    f_fft = rfftfreq(n_steps + 1, 1000 / Δt)

    max_f = maximum(abs ∘ f_analytic, k_fft)

    f_ind_range = filter(i -> abs(f_fft[i]) < max_f * 1.5, eachindex(f_fft))
    f_fft = @view f_fft[f_ind_range]
    Sx_fft = @view Sx_fft[:, f_ind_range]

    max_amp = maximum(norm, Sx_fft)

    f, ax, _ = heatmap(k_fft, f_fft, norm.(Sx_fft);
        colorrange=(max_amp * cutoff, max_amp / gain))
    lines!(ax, k_fft, f_analytic; color=:red)

    limits!(ax, -maximum(k_fft), maximum(k_fft), 0, max_f * 1.1)

    ax.xlabel[] = "k"
    ax.ylabel[] = "f [1/ps]"

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

function test_magnetization(n)
    S = randn(n, n, n, 3)
    normalize_spin!(S)

    state = init_state_par(S)
    params = setup_params(
        20.0, 3.0, 0.0, (@SVector [0.0, 0.0, 0.0]), 1.0, 0.1
    )

    M_hist = @time simulate_magnetization!(state, params, 100000)
    lines(M_hist)
end

function test_demagnetization(n)
    S = zeros(n, n, n, 3)
    S[:, :, :, 3] .= 1.0

    state = init_state_par(S)
    params = setup_params(
        10.0, 3.0, 10.0, (@SVector [0.0, 0.0, 0.0]), 1.0, 0.1
    )

    n_steps = 1000
    M_hist = @time simulate_magnetization!(state, params, n_steps)
    f, _, _ = lines(M_hist)
    ax2 = Axis(f[2, 1])

    data_points = @view M_hist[n_steps÷100:end]

    pts = 0:100:length(data_points)
    lines!(ax2, pts, autocor(data_points, pts))

    M_mean = mean(data_points)
    M_stdd = std(data_points)

    @show M_mean M_stdd

    uncorr_spacing = n_steps ÷ 100
    data_points_uncorr = @view data_points[1:uncorr_spacing:end]
    M_mean = mean(data_points_uncorr)
    M_stdd = std(data_points_uncorr)

    @show M_mean M_stdd

    scatter!(ax2, [uncorr_spacing], autocor(data_points, [uncorr_spacing]);
        color=:red)

    f
end

function make_phase_diagram(n)
    T_range = range(0.0, 20.0, 200)

    M_means = [0.0 for _ in T_range]
    M_stdds = [0.0 for _ in T_range]

    normal_steps = 10000
    critical_steps = 100 * normal_steps
    critical_temp = 16.0

    steps_func(kT) = ceil(Int, normal_steps +
                               (critical_steps - normal_steps) *
                               exp(-(kT - critical_temp)^2))

    total_steps = sum(steps_func, T_range)

    progress_meter = Channel{Int}(spawn=true) do ch
        done_steps = 0
        while done_steps < total_steps
            steps = take!(ch)
            done_steps += steps
            println("$done_steps / $total_steps = ",
                round(100 * done_steps / total_steps, digits=2), "%")
        end
    end

    Threads.@threads for i in eachindex(T_range)
        kT = T_range[i]

        params = setup_params(
            10.0, 3.0, kT, (@SVector [0.0, 0.0, 0.0]), 1.0, 0.1
        )

        steps = steps_func(kT)
        # @show kT, steps

        M_mean, M_stdd = simulate_demagnetization(
            n, steps, 0.1, params
        )

        M_means[i] = M_mean
        M_stdds[i] = M_stdd

        put!(progress_meter, steps)
    end

    f, ax, _ = band(T_range, M_means - M_stdds, M_means + M_stdds)
    lines!(ax, T_range, M_means)

    f
end
