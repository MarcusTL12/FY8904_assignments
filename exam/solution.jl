using LinearAlgebra
import Plots
using Printf
using Statistics
using Distributed

include("visualization.jl")
include("lattice2d.jl")
include("lattice3d.jl")
include("mc2d.jl")
include("mc3d.jl")
include("task1.jl")
include("task2.jl")

# This function initializes the Distributed.jl package with the number of
# workers requested. If this is not run, the phase diagrams will run single
# threaded, which can take a while.
function setup_distributed(n_workers)
    w = workers()
    if w == [1]
        w = addprocs(1)
    end
    if length(w) < n_workers
        addprocs(n_workers - length(w))
        @everywhere include("solution.jl")
    elseif length(w) > n_workers
        Δ = length(w) - n_workers
        rmprocs(w[end-(Δ-1):end])
    end
end
