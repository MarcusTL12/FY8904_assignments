using Plots
using LsqFit
using Measurements

function simulate_walk(dist, n)
    xs = sizehint!(Float64[], n)

    x = 0.0

    for _ in 1:n
        dx = dist()
        x += dx
        push!(xs, x)
    end

    xs
end

function simulate_first_return(dist, max_steps)
    x = 0.0

    for t in 0:max_steps
        dx = dist()
        nx = x + dx
        if x * nx < 0
            return t
        end
        x = nx
    end

    max_steps
end

function simulate_n_returns(dist, max_steps, n)
    ts = zeros(Int, n)

    Threads.@threads for i in 1:n
        ts[i] = simulate_first_return(dist, max_steps)
    end

    ts
end

function make_counts(ts, max_steps, counts=zeros(Int, max_steps))
    for t in ts
        if t <= length(counts)
            counts[t] += 1
        end
    end

    counts
end

function make_prob_dist(counts, n_simulations)
    probs = zeros(Float64, length(counts))

    for (i, c) in enumerate(counts)
        probs[i] = c / n_simulations
    end

    probs
end

function run_a()
    max_steps = 1000
    n_simulations = 100_000_000

    ts = simulate_n_returns(randn, max_steps + 1, n_simulations)
    cs = make_counts(ts, max_steps)
    ps = make_prob_dist(cs, n_simulations)

    plot((@view ps[1:100]); leg=false)
end

function find_alpha(ps, p0=[1.0, 1.0])
    model(t, (c, a)) = @. c * t^-a

    fit = curve_fit(
        model, collect(eachindex(ps)), ps, p0;
        autodiff=:forwarddiff
    )

    fit.param, stderror(fit)
end

function run_b(max_steps, n_simulations)
    ts = simulate_n_returns(randn, max_steps + 1, n_simulations)
    cs = make_counts(ts, max_steps)
    ps = make_prob_dist(cs, n_simulations)

    p, e = find_alpha(ps)

    println("α = ", p[2] ± e[2])

    p, cs
end

function improve_b((p, cs), n_simulations)
    max_steps = length(cs)

    ts = simulate_n_returns(randn, max_steps + 1, n_simulations)
    make_counts(ts, max_steps, cs)
    ps = make_prob_dist(cs, n_simulations)

    p, e = find_alpha(ps, p)

    println("α = ", p[2] ± e[2])

    p, cs
end
