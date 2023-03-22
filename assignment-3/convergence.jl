
function converge_n(l, ns, nev)
    es = Union{Float64,Missing}[]

    for n in ns
        lattice = make_lattice(l, n)
        A = make_linmap(lattice, l, n)
        @show A.N
        @time if A.N <= nev
            B = make_reduced_matrix(lattice, get_h(l, n))
            e = eigvals(Symmetric(Matrix(B)))
            append!(es, e)
            append!(es, missing for _ in (length(e)+1):nev)
        else
            e, _ = eigs(A; nev=nev, which=:SR, maxiter=10000)
            append!(es, e)
        end
    end

    for (i, e) in enumerate(es)
        es[i] = √e
    end

    es = reshape(es, (nev, length(ns)))
end

function converge_ln(ls, ns, nev)
    es = Union{Float64,Missing}[]

    for (l, n) in zip(ls, ns)
        lattice = make_lattice(l, n)
        A = make_linmap(lattice, l, n)
        @show A.N
        @time if A.N <= nev
            B = make_reduced_matrix(lattice, get_h(l, n))
            e = eigvals(Symmetric(Matrix(B)))
            append!(es, e)
            append!(es, missing for _ in (length(e)+1):nev)
        else
            e, _ = eigs(A; nev=nev, which=:SR, maxiter=10000)
            append!(es, e)
        end
    end

    for (i, e) in enumerate(es)
        es[i] = √e
    end

    es = reshape(es, (nev, length(ns)))
end

function plot_convergence(es)
    Plots.plot(; leg=false)

    for r in eachrow(es)
        Plots.plot!(r)
    end

    Plots.plot!()
end
