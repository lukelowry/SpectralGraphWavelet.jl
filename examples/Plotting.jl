using SparseArrays, LinearAlgebra, SuiteSparse, Plots

include("../src/GraphFilters.jl")
include("../src/GraphIO.jl")
include("../src/GraphPlot.jl")

base_dir = joinpath(@__DIR__, "..", "..", "python", "sgwt", "library", "data")
path_L = joinpath(base_dir, "DELAY", "TEXAS_DELAY.mat")
path_S = joinpath(base_dir, "SIGNALS", "texas_coords.mat")

A = load_laplacian(path_L)
S = load_signal(path_S)

scale = 0.1
conv = DyConvolve(A, [1/scale]) # pass poles (1/s)
b = impulse(A, 1201) + impulse(A, 601)

results_bp = bandpass(conv, b)[1]

f_vals = results_bp 
plot_signal(f_vals, S)

println("Press Enter to exit...")
readline()