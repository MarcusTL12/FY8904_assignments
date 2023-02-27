
function space_averaged_magnetization(S)
    mean(@view S[:, :, :, 3])
end

function simulate_demagnetization(n, n_step, eq_frac, params::SimParams)
    S = zeros(n, n, n, 3)
    S[:, :, :, 3] .= 1.0

    state = init_state_par(S)

    M_hist = simulate_magnetization!(state, params, n_step)

    eq_steps = ceil(Int, eq_frac * n_step)

    data_points = @view M_hist[eq_steps:end]

    mean(data_points), std(data_points)
end
