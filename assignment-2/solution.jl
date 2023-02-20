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
