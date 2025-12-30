include("../src/GraphFilters.jl")
include("../src/GraphIO.jl")
include("../src/GraphPlot.jl")

base_dir = joinpath(@__DIR__, "..", "..", "python", "sgwt", "library", "data")
path_L = joinpath(base_dir, "DELAY", "USA_DELAY.mat")
path_S = joinpath(base_dir, "SIGNALS", "usa_coords.mat")

A = load_laplacian(path_L)
S = load_signal(path_S)


scales = 10 .^ range(-5, 2, length=20)
poles = 1 ./scales

conv = DyConvolve(A, poles)
b = impulse(A, 1201) 


println("Warming up JIT compiler...")

# 1. Warmup: Run once to trigger compilation (result ignored)
bandpass(conv, b)

println("Benchmarking...")

@time begin
    for i in 1:2
        global results_bp = bandpass(conv, b)
    end
end

println("Done!")