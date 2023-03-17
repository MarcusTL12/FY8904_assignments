using LinearAlgebra

function make_mat(n)
    A = diagm(0 => randn(n), 1 => randn(n - 1), n ÷ 2 => randn(n - n ÷ 2))
    A = Symmetric(A)
end

function make_symmat(n)
    SymDiags(0 => randn(n), 1 => randn(n - 1), n ÷ 2 => randn(n - n ÷ 2))
end

# Symmetric n×n matrix with a few diagonals
struct SymDiags{T} <: AbstractMatrix{T}
    n::Int
    diagonals::Vector{Tuple{Int,Vector{T}}}

    function SymDiags(kv::Pair{Int,Vector{T}}...) where {T}
        diagonals = sort!(collect((k, v) for (k, v) in kv))

        (f_i, f_v), rest = Iterators.peel(kv)
        n = length(f_v) - f_i

        for (i, v) in rest
            @assert length(v) + i == n
        end

        new{T}(n, diagonals)
    end
end

Base.size(A::SymDiags) = (A.n, A.n)

function Base.getindex(A::SymDiags, i::Int)
    A[fldmod1(i, A.n)...]
end

function Base.getindex(A::SymDiags{T}, i::Int, j::Int) where {T}
    if i > j
        i, j = j, i
    end

    diag_n = j - i

    diag_i = findfirst(x -> x[1] == diag_n, A.diagonals)

    if !isnothing(diag_i)
        A.diagonals[diag_i][2][i]
    else
        zero(T)
    end
end
