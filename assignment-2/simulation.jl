using Random
using StaticArrays

mutable struct SimState
    S::Array{Float64,4}
    Sp::Array{Float64,4}
    ∂S::Array{Float64,4}
    ∂Sp::Array{Float64,4}
    Γ1::Array{Float64,4}
    Γ2::Array{Float64,4}
end

struct SimParams
    J::Float64
    dz::Float64
    B::SVector{3,Float64}
    α::Float64
    γ::Float64
    μ::Float64
    kT::Float64
    Δt::Float64
end

function init_state(S)
    SimState(S, similar(S), similar(S), similar(S), randn(size(S)), similar(S))
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
        params.γ, params.μ, params.kT, params.Δt, Γ)
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

function do_heun_step!(state::SimState, params::SimParams)
    Γ_task = Threads.@spawn randn!(state.Γ2)

    # f(tn, yn)
    compute_∂S!(state.∂S, state.S, state.Γ1, params)

    # ypn+1
    @inline do_normalized_euler_step!(state.Sp, state.S, state.∂S, params.Δt)

    wait(Γ_task)

    # f(tn+1, ypn+1)
    compute_∂S!(state.∂Sp, state.Sp, state.Γ2, params)

    # 0.5 Δt (f(tn, yn) + f(tn+1, ypn+1))
    c = 0.5 * params.Δt
    @turbo for i in eachindex(state.S)
        state.∂S[i] = c * (state.∂S[i] + state.∂Sp[i])
    end

    # yn+1
    @inline do_normalized_euler_step!(state.S, state.S, state.∂S, 1.0)

    state.Γ1, state.Γ2 = state.Γ2, state.Γ1
end
