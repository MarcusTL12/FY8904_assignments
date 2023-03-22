
function plot_idos(e, l)
    ω = sqrt.(e)
    N = [i - 1 for i in eachindex(e)]

    Plots.plot(ω, N)

    analytic_N = [1 / 4π * ω^2 - 2^l / 4π * ω for ω in ω]

    Plots.plot!(ω, analytic_N)
end

function get_Δidos(e)
    ω = sqrt.(e)
    N = [i - 1 for i in eachindex(e)]

    @. ω^2 / 4π - N
end

function plot_Δidos(e)
    ω = sqrt.(e)
    ΔN = get_Δidos(e)

    Plots.plot(ω, ΔN)
end
