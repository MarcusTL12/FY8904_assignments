Here is a quick guide to the contents of this zip file.

The main entry point of the code is the solution.jl file. To run the code
one would open a julia terminal and include this file:

julia> include("solution.jl")

then you can run the specific function you want, for example:

julia> run_2_1_7a()

The core logic of the 2d simulations is in the lattice2d.jl and mc2d.jl files.
A lot of the functions in these files also happened to work for 3d. The
3d specific versions of the functions that did not work automatically are in
the lattice3d.jl and mc3d.jl files.

In the task1.jl and task2.jl files, there are all the functions you would
call to perform the simulations for a given task, for example run_2_2_3()
performs the simulations to plot the phase diagrams for a 3d polymer.

Interactive plots:
There are two interactive plots included; chain_100.html and chain_200.html.
These can be opened in a browser to view.
