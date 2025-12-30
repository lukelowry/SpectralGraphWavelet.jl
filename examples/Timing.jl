using SpectralGraphWavelet

script_dir = dirname(@__FILE__)

path_L = joinpath(script_dir, "USA_DELAY.mat")
path_S = joinpath(script_dir, "usa_coords.mat")

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