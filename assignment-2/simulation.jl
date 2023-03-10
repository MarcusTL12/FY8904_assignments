
abstract type SimState end

mutable struct SimStateSerial <: SimState
    S::Array{Float64,4}
    Sp::Array{Float64,4}
    ∂S::Array{Float64,4}
    ∂Sp::Array{Float64,4}
    Γ1::Array{Float64,4}
    Γ2::Array{Float64,4}
end

mutable struct SimStateParallel <: SimState
    S::Array{Float64,4}
    Sp::Array{Float64,4}
    ∂S::Array{Float64,4}
    ∂Sp::Array{Float64,4}
    Γ::Array{Float64,4}

    Γ_send::Channel{Array{Float64,4}}
    Γ_recv::Channel{Array{Float64,4}}
end

struct SimParams
    J::Float64
    dz::Float64
    B::SVector{3,Float64}
    α::Float64
    δ::Float64
    kT::Float64
    Δt::Float64
end

# This assumes that parameters are given in the following units:
# J: meV
# dz: meV
# kT: meV
# B: meV x 3
# Δt: fs
# α: unitless
function setup_params(J, dz, kT, B, Δt, α)
    # defining δ = γ/μ
    # Here is δ in units [1 / (fs meV)]
    δ = 0.0027586206896551726

    SimParams(J, dz, B, α, δ, kT, Δt)
end

function init_state(S)
    SimStateSerial(
        S, similar(S), similar(S), similar(S), randn(size(S)), similar(S)
    )
end

function init_state_par(S)
    Γ_recv = Channel{Array{Float64,4}}(1)

    Γ_send = Channel{Array{Float64,4}}(1; spawn=true) do ch
        while true
            Γ = take!(ch)
            randn!(Γ)
            put!(Γ_recv, Γ)
        end
    end

    Γ2 = similar(S)
    put!(Γ_send, Γ2)

    state = SimStateParallel(
        S, similar(S), similar(S), similar(S), randn(size(S)), Γ_send, Γ_recv
    )

    finalizer(state) do state
        close(state.Γ_send)
        close(state.Γ_recv)
    end

    state
end

function normalize_spin!(S)
    @inbounds @fastmath begin
        Sxs = @view S[:, :, :, 1]
        Sys = @view S[:, :, :, 2]
        Szs = @view S[:, :, :, 3]

        @turbo for i in eachindex(Sxs)
            Sx = Sxs[i]
            Sy = Sys[i]
            Sz = Szs[i]

            l_inv = 1 / √(Sx^2 + Sy^2 + Sz^2)

            Sxs[i] = Sx * l_inv
            Sys[i] = Sy * l_inv
            Szs[i] = Sz * l_inv
        end
    end
end

function compute_∂S!(∂S, S, Γ, params::SimParams)
    @inline compute_∂S!(∂S, S, params.J, params.dz, params.B, params.α,
        params.δ, params.kT, params.Δt, Γ)
end

function do_normalized_euler_step!(S2, S1, ∂S, c)
    @inbounds @fastmath begin
        S2xs = @view S2[:, :, :, 1]
        S2ys = @view S2[:, :, :, 2]
        S2zs = @view S2[:, :, :, 3]

        S1xs = @view S1[:, :, :, 1]
        S1ys = @view S1[:, :, :, 2]
        S1zs = @view S1[:, :, :, 3]

        ∂Sxs = @view ∂S[:, :, :, 1]
        ∂Sys = @view ∂S[:, :, :, 2]
        ∂Szs = @view ∂S[:, :, :, 3]

        @turbo for i in eachindex(S1xs)
            Sx = S1xs[i] + c * ∂Sxs[i]
            Sy = S1ys[i] + c * ∂Sys[i]
            Sz = S1zs[i] + c * ∂Szs[i]

            l_inv = 1 / √(Sx^2 + Sy^2 + Sz^2)

            S2xs[i] = Sx * l_inv
            S2ys[i] = Sy * l_inv
            S2zs[i] = Sz * l_inv
        end
    end
end

function do_heun_step!(state::SimStateSerial, params::SimParams)
    # f(tn, yn)
    compute_∂S!(state.∂S, state.S, state.Γ1, params)

    # ypn+1
    @inline do_normalized_euler_step!(state.Sp, state.S, state.∂S, params.Δt)

    # Get new random numbers
    randn!(state.Γ2)

    # f(tn+1, ypn+1)
    compute_∂S!(state.∂Sp, state.Sp, state.Γ2, params)

    # 0.5 Δt (f(tn, yn) + f(tn+1, ypn+1))
    c = 0.5 * params.Δt
    @turbo for i in eachindex(state.S)
        state.∂S[i] = c * (state.∂S[i] + state.∂Sp[i])
    end

    # yn+1
    @inline do_normalized_euler_step!(state.S, state.S, state.∂S, 1.0)

    # Swap random numbers
    Γ_tmp = state.Γ1
    state.Γ1 = state.Γ2
    state.Γ2 = Γ_tmp
end

function do_heun_step!(state::SimStateParallel, params::SimParams)
    # f(tn, yn)
    compute_∂S!(state.∂S, state.S, state.Γ, params)
    put!(state.Γ_send, state.Γ) # Send off Γ to be randomized

    # ypn+1
    @inline do_normalized_euler_step!(state.Sp, state.S, state.∂S, params.Δt)

    state.Γ = take!(state.Γ_recv) # Take randomized Γ

    # f(tn+1, ypn+1)
    compute_∂S!(state.∂Sp, state.Sp, state.Γ, params)

    # 0.5 Δt (f(tn, yn) + f(tn+1, ypn+1))
    c = 0.5 * params.Δt
    @turbo for i in eachindex(state.S)
        state.∂S[i] = c * (state.∂S[i] + state.∂Sp[i])
    end

    # yn+1
    @inline do_normalized_euler_step!(state.S, state.S, state.∂S, 1.0)
end

function simulate!(state::SimState, params::SimParams,
    n_frames, n_steps_per_frame, S_hist=Float64[], add_first=true)

    S_hist = reshape(S_hist, length(S_hist))

    if add_first
        append!(S_hist, state.S)
    end

    for _ in 1:n_frames
        for _ in 1:n_steps_per_frame
            do_heun_step!(state, params)
        end

        append!(S_hist, state.S)
    end

    nx, ny, nz, _ = size(state.S)
    n = nx * ny * nz

    reshape(S_hist, n, 3, length(S_hist) ÷ 3n)
end

function simulate_magnetization!(state::SimState, params::SimParams, n_steps)
    M_hist = [space_averaged_magnetization(state.S)]

    for _ in 1:n_steps
        do_heun_step!(state, params)
        push!(M_hist, space_averaged_magnetization(state.S))
    end

    M_hist
end

function simulate_2d_surface!(state::SimState, params::SimParams,
    n_frames, n_steps_per_frame)

    nx, ny, nz, _ = size(state.S)
    @assert nz == 1

    anim_threads = max(1, Threads.nthreads() - 2)
    dir = "tmp_frames"

    rm(dir; recursive=true, force=true)
    mkdir(dir)

    buffer_ch = Channel{Array{Float64,4}}(anim_threads + 2)

    for _ in 1:anim_threads+1
        put!(buffer_ch, similar(state.S))
    end

    workers = [
        Channel{Tuple{Int,Array{Float64,4}}}(; spawn=false) do work
            while true
                i, S = take!(work)
                Sr = reshape(S, nx, ny, 3)
                img = visualize_spin_surface(Sr)
                save("$dir/$i.png", img)
                put!(buffer_ch, S)
            end
        end for _ in 1:anim_threads
    ]

    w = 1

    for i in 1:n_frames
        for _ in 1:n_steps_per_frame
            do_heun_step!(state, params)
        end

        buffer = take!(buffer_ch)
        copy!(buffer, state.S)
        put!(workers[w], (i, buffer))
        w += 1
        if w > anim_threads
            w = 1
        end
    end

    # close(worker)
    for worker in workers
        close(worker)
    end

    make_mp4(dir, "tmp_anim")
end
