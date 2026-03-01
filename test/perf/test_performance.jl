using SpectralGraphWavelet
using SparseArrays

function path_laplacian_perf(n::Int)
    d0 = fill(2.0, n)
    d0[1] = 1.0
    d0[end] = 1.0
    return spdiagm(-1 => fill(-1.0, n - 1), 0 => d0, 1 => fill(-1.0, n - 1))
end

function run_perf_suite()
    println("Running SGWT performance smoke suite...")
    for n in (100, 500, 2000)
        L = path_laplacian_perf(n)
        x = randn(n)
        scales = [1.0]
        conv = Convolve(L)
        t_lp = @elapsed lowpass(conv, x, scales)
        t_bp = @elapsed bandpass(conv, x, scales)
        t_hp = @elapsed highpass(conv, x, scales)
        println("n=$n  lowpass=$(round(t_lp, digits=4))s  bandpass=$(round(t_bp, digits=4))s  highpass=$(round(t_hp, digits=4))s")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_perf_suite()
end
