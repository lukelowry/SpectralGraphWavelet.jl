using SparseArrays

mutable struct ChebyConvolve
    n_vertices::Int
    L::SparseMatrixCSC{Float64, Int}
    spectrum_bound::Float64
    recurrence_cache::Dict{Tuple{Float64, Float64}, SparseMatrixCSC{Float64, Int}}
end

function ChebyConvolve(L::SparseMatrixCSC)
    A = SparseMatrixCSC{Float64, Int}(L)
    ub = estimate_spectral_bound(A)
    return ChebyConvolve(size(A, 1), A, ub, Dict{Tuple{Float64, Float64}, SparseMatrixCSC{Float64, Int}}())
end
ChebyConvolve(L::AbstractSparseMatrix) = ChebyConvolve(SparseMatrixCSC{Float64, Int}(L))

function _get_recurrence_matrix(conv::ChebyConvolve, spectrum_bound::Float64, min_lambda::Float64=0.0)
    key = (spectrum_bound, min_lambda)
    if haskey(conv.recurrence_cache, key)
        return conv.recurrence_cache[key]
    end
    range_ = spectrum_bound - min_lambda
    alpha = 2.0 / range_
    beta = -(spectrum_bound + min_lambda) / range_
    M = alpha .* conv.L .+ beta .* spdiagm(0 => ones(Float64, conv.n_vertices))
    conv.recurrence_cache[key] = M
    return M
end

function _accumulate!(W::Array{Float64, 3}, Z::AbstractMatrix{Float64}, c::AbstractVector{Float64})
    nd = length(c)
    if nd <= 4
        for d in 1:nd
            c[d] == 0 && continue
            W[:, :, d] .+= Z .* c[d]
        end
    else
        for d in 1:nd
            W[:, :, d] .+= Z .* c[d]
        end
    end
    return W
end

function _convolve_real(conv::ChebyConvolve, B::AbstractMatrix{Float64}, C::ChebyKernel)
    n_vertex, n_signals = size(B)
    n_order, n_dim = size(C.C)
    if n_order == 0 || n_dim == 0
        return zeros(Float64, n_vertex, n_signals, n_dim)
    end

    W = zeros(Float64, n_vertex, n_signals, n_dim)
    coeffs = Matrix{Float64}(C.C)

    Tkm2 = copy(B)
    _accumulate!(W, Tkm2, vec(coeffs[1, :]))
    if n_order == 1
        return W
    end

    M = _get_recurrence_matrix(conv, C.spectrum_bound, C.min_lambda)
    Tkm1 = M * Tkm2
    _accumulate!(W, Tkm1, vec(coeffs[2, :]))
    for k in 3:n_order
        Tk = 2.0 .* (M * Tkm1) .- Tkm2
        Tkm2 = Tkm1
        Tkm1 = Tk
        _accumulate!(W, Tkm1, vec(coeffs[k, :]))
    end
    return W
end

function convolve(conv::ChebyConvolve, B, C::ChebyKernel)
    if any(!isreal, B)
        was_1d = ndims(B) == 1
        Bmat = was_1d ? reshape(B, :, 1) : B
        Wre = _convolve_real(conv, Float64.(real(Bmat)), C)
        Wim = _convolve_real(conv, Float64.(imag(Bmat)), C)
        W = Wre + im * Wim
        if was_1d
            return size(W, 2) == 1 ? dropdims(W; dims=2) : W
        end
        return W
    end

    was_1d = ndims(B) == 1
    Bmat = was_1d ? reshape(Float64.(B), :, 1) : Float64.(B)
    W = _convolve_real(conv, Bmat, C)
    if was_1d
        return size(W, 2) == 1 ? dropdims(W; dims=2) : W
    end
    return W
end

function convolve_multi(conv::ChebyConvolve, B, kernels::AbstractVector)
    isempty(kernels) && return Any[]
    return [convolve(conv, B, k) for k in kernels]
end
