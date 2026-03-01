using LinearAlgebra
using Random
using SparseArrays

function path_laplacian(n::Int)
    d0 = fill(2.0, n)
    d0[1] = 1.0
    d0[end] = 1.0
    return spdiagm(-1 => fill(-1.0, n - 1), 0 => d0, 1 => fill(-1.0, n - 1))
end

function random_signal(L::SparseMatrixCSC; ncols::Int=5, seed::Int=42)
    Random.seed!(seed)
    return randn(size(L, 1), ncols)
end

function first_available_laplacian()
    names = [k for k in keys(SpectralGraphWavelet._ensure_registry()) if startswith(k, "DELAY_")]
    isempty(names) && error("No bundled DELAY_* resources discovered.")
    return resource(first(sort(names)))
end
