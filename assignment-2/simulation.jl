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

function compute_∂S!(∂S, S, Γ, params::SimParams)
    compute_∂S!(∂S, S, params.J, params.dz, params.B, params.α, params.γ,
        params.μ, params.kT, params.Δt, Γ)
end

function do_heun_step!(state::SimState, params::SimParams)
    # f(tn, yn)
    compute_∂S!(state.∂S, state.S, state.Γ1, params)

    # ypn+1
    @inbounds @simd for i in eachindex(state.S)
        state.Sp[i] = state.S[i] + state.∂S[i] * params.Δt
    end

    randn!(state.Γ2)

    # f(tn+1, ypn+1)
    compute_∂S!(state.∂Sp, state.Sp, state.Γ2, params)

    # 0.5 Δt (f(tn, yn) + f(tn+1, ypn+1))
    c = 0.5 * params.Δt
    for i in eachindex(state.S)
        state.∂S[i] = c * (state.∂S[i] + state.∂Sp[i])
    end

    # yn+1
    @inbounds @simd for i in eachindex(state.S)
        state.S[i] = state.S[i] + state.∂S[i]
    end

    state.Γ1, state.Γ2 = state.Γ2, state.Γ1
end
