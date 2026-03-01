using LinearAlgebra
using Statistics

struct ChebyModel
    C::Matrix{Float64}
    spectrum_bound::Float64
    min_lambda::Float64
end

function (m::ChebyModel)(x)
    xv = Float64.(vec(x))
    x_scaled = @. (2.0 * (xv - m.min_lambda) / (m.spectrum_bound - m.min_lambda)) - 1.0
    y = _chebval_matrix(x_scaled, m.C)
    return size(y, 2) == 1 ? vec(y) : y
end

function _ensure_2d(arr)
    a = Float64.(arr)
    if ndims(a) == 1
        return reshape(a, :, 1)
    end
    return a
end

function _truncate(coeffs::AbstractMatrix, rtol::Real)
    thr = rtol * maximum(abs.(coeffs))
    rowmax = vec(maximum(abs.(coeffs); dims=2))
    significant = findall(>(thr), rowmax)
    lastidx = isempty(significant) ? 1 : (last(significant))
    return Matrix{Float64}(coeffs[1:lastidx, :])
end

function _compute_coefficients(f, order::Int, spectrum_bound::Float64, min_lambda::Float64, n_samples, sampling::String)
    lambda_range = spectrum_bound - min_lambda
    lambda_mid = (spectrum_bound + min_lambda) / 2.0

    if sampling == "chebyshev"
        n = order + 1
        k = collect(0:order)
        x_cheb = cos.(pi .* k ./ order)
        sample_x = lambda_mid .+ (lambda_range / 2.0) .* x_cheb
        fvals = _ensure_2d(f(sample_x))
        coeffs = zeros(Float64, n, size(fvals, 2))
        w = ones(Float64, n)
        w[1] = 0.5
        w[end] = 0.5
        for j in 0:order
            Tj = cos.(j * pi .* k ./ order)
            scale = (j == 0 || j == order) ? (1.0 / order) : (2.0 / order)
            coeffs[j + 1, :] = scale .* vec(sum((w .* Tj) .* fvals; dims=1))
        end
        return coeffs
    end

    ns = n_samples === nothing ? max(4 * (order + 1), 1000) : n_samples
    t = range(0.0, 1.0, length=ns)
    sample_x = if sampling == "quadratic"
        min_lambda .+ lambda_range .* (t .^ 2)
    elseif sampling == "logarithmic"
        eps = max(min_lambda * 0.001, 1e-10)
        exp.(log(min_lambda + eps) .+ t .* log(spectrum_bound / (min_lambda + eps)))
    else
        min_lambda .+ lambda_range .* t
    end
    x_scaled = @. 2.0 * (sample_x - min_lambda) / lambda_range - 1.0
    fvals = _ensure_2d(f(sample_x))

    # Weighted least-squares fit using Vandermonde of Chebyshev basis.
    weights = @. 1.0 / sqrt(1.0 - min(x_scaled^2, 0.9999))
    T = zeros(Float64, length(x_scaled), order + 1)
    T[:, 1] .= 1.0
    if order >= 1
        T[:, 2] .= x_scaled
    end
    for k in 3:(order + 1)
        T[:, k] .= 2.0 .* x_scaled .* T[:, k - 1] .- T[:, k - 2]
    end
    W = Diagonal(weights)
    coeffs = (W * T) \ (W * fvals)
    return coeffs
end

function _adaptive_fit(f, start_order::Int, spectrum_bound::Float64, min_lambda::Float64, sampling::String, rtol::Float64, max_order::Int, target_error::Float64)
    test_x = collect(range(min_lambda, spectrum_bound, length=1000))
    f_exact = _ensure_2d(f(test_x))
    order = max(start_order, 8)
    coeffs = zeros(1, 1)
    while true
        coeffs = _ensure_2d(_compute_coefficients(f, order, spectrum_bound, min_lambda, nothing, sampling))
        x_scaled = @. 2.0 * (test_x - min_lambda) / (spectrum_bound - min_lambda) - 1.0
        f_approx = _chebval_matrix(x_scaled, coeffs)
        rel = maximum(abs.(f_exact .- f_approx) ./ max.(abs.(f_exact), 1e-15))
        if rel <= target_error || order >= max_order
            break
        end
        order = min(Int(floor(order * 1.5) + 1), max_order)
    end
    return ChebyModel(_truncate(coeffs, rtol), spectrum_bound, min_lambda)
end

function fit(
    ::Type{ChebyModel},
    f::Function,
    order::Int,
    spectrum_bound::Float64;
    n_samples=nothing,
    sampling::String="chebyshev",
    min_lambda::Float64=0.0,
    rtol::Float64=1e-12,
    adaptive::Bool=false,
    max_order::Int=500,
    target_error::Float64=1e-10,
)
    order < 1 && throw(ArgumentError("Order must be >= 1"))
    if adaptive
        return _adaptive_fit(f, order, spectrum_bound, min_lambda, sampling, rtol, max_order, target_error)
    end
    coeffs = _ensure_2d(_compute_coefficients(f, order, spectrum_bound, min_lambda, n_samples, sampling))
    coeffs = _truncate(coeffs, rtol)
    return ChebyModel(coeffs, spectrum_bound, min_lambda)
end

function kernel(::Type{ChebyModel}, L::SparseMatrixCSC, f::Function, order::Int; kwargs...)
    ub = estimate_spectral_bound(L)
    model = fit(ChebyModel, f, order, ub; kwargs...)
    return ChebyKernel(C=model.C, spectrum_bound=model.spectrum_bound, min_lambda=model.min_lambda)
end

function kernel_on_graph(::Type{ChebyModel}, L::SparseMatrixCSC, f::Function, order::Int; kwargs...)
    return kernel(ChebyModel, L, f, order; kwargs...)
end

function Base.getproperty(::Type{ChebyModel}, s::Symbol)
    if s === :fit
        return (args...; kwargs...) -> fit(ChebyModel, args...; kwargs...)
    elseif s === :kernel
        return (args...; kwargs...) -> kernel(ChebyModel, args...; kwargs...)
    elseif s === :kernel_on_graph
        return (args...; kwargs...) -> kernel_on_graph(ChebyModel, args...; kwargs...)
    end
    return getfield(ChebyModel, s)
end
