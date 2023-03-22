using Plots
using DataStructures
using Arpack
using LinearMaps
using SparseArrays
using KrylovKit
using ArnoldiMethod

include("fractal.jl")
include("lattice.jl")
include("5point-stencil.jl")
include("analysis.jl")
include("anim.jl")
