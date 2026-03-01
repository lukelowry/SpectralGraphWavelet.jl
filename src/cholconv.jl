using LinearAlgebra
using SparseArrays
using SuiteSparse

function _squeeze_result(x)
    if ndims(x) == 2 && size(x, 2) == 1
        return vec(x)
    elseif ndims(x) == 3 && size(x, 2) == 1
        return dropdims(x; dims=2)
    end
    return x
end

function _process_signal(func::Function, B_in; scales=nothing)
    was_1d = ndims(B_in) == 1
    B = was_1d ? reshape(B_in, :, 1) : B_in

    scalar_scale = scales !== nothing && scales isa Real
    scales_vec = if scales === nothing
        nothing
    elseif scalar_scale
        [Float64(scales)]
    else
        Float64.(collect(scales))
    end

    if any(!isreal, B)
        Breal = Float64.(real(B))
        Bimag = Float64.(imag(B))
        real_part = scales_vec === nothing ? func(Breal) : func(Breal, scales_vec)
        imag_part = scales_vec === nothing ? func(Bimag) : func(Bimag, scales_vec)
        result = if real_part isa AbstractVector
            [r + im * i for (r, i) in zip(real_part, imag_part)]
        else
            real_part + im * imag_part
        end
    else
        Bf = Float64.(B)
        result = scales_vec === nothing ? func(Bf) : func(Bf, scales_vec)
    end

    if was_1d
        if result isa AbstractVector
            result = [_squeeze_result(r) for r in result]
        else
            result = _squeeze_result(result)
        end
    end

    if scalar_scale && result isa AbstractVector
        result = result[1]
    end
    return result
end

mutable struct Convolve
    n_vertices::Int
    L::SparseMatrixCSC{Float64, Int}
    factor_cache::Dict{Float64, SuiteSparse.CHOLMOD.Factor{Float64}}
end

function Convolve(L::SparseMatrixCSC)
    A = SparseMatrixCSC{Float64, Int}(L)
    return Convolve(size(A, 1), A, Dict{Float64, SuiteSparse.CHOLMOD.Factor{Float64}}())
end
Convolve(L::AbstractSparseMatrix) = Convolve(SparseMatrixCSC{Float64, Int}(L))

function _factor!(conv::Convolve, q::Float64; force::Bool=true)
    if !force && haskey(conv.factor_cache, q)
        return conv.factor_cache[q]
    end
    F = ldlt(Symmetric(conv.L); shift=q)
    conv.factor_cache[q] = F
    return F
end

function _as_vfkernel(K)
    if K isa VFKernel
        return K
    elseif K isa AbstractDict
        return VFKernel_from_dict(K)
    end
    throw(ArgumentError("Kernel K must be a VFKernel object or a compatible dictionary."))
end

function _kernel_arrays(K::VFKernel)
    if K.R === nothing || K.Q === nothing
        throw(ArgumentError("Kernel K must contain residues (R) and poles (Q)."))
    end
    R = K.R isa AbstractVector ? reshape(Float64.(K.R), :, 1) : Float64.(Matrix(K.R))
    Q = Float64.(vec(K.Q))
    D = if K.D === nothing
        Float64[]
    elseif K.D isa Number
        Float64[K.D]
    else
        Float64.(vec(K.D))
    end
    return R, Q, D
end

function _accum_residue!(W::Array{Float64, 3}, Z::AbstractMatrix{Float64}, r::AbstractVector{Float64})
    for d in eachindex(r)
        W[:, :, d] .+= Z .* r[d]
    end
    return W
end

function _convolve_impl(conv::Convolve, B::AbstractMatrix{Float64}, K)
    kern = _as_vfkernel(K)
    R, Q, D = _kernel_arrays(kern)
    nDim = size(R, 2)
    W = zeros(Float64, size(B, 1), size(B, 2), nDim)
    if !isempty(D)
        for d in eachindex(D)
            W[:, :, d] .+= B .* D[d]
        end
    end
    for (i, q) in enumerate(Q)
        F = _factor!(conv, q; force=true)
        Z = F \ B
        _accum_residue!(W, Z, vec(R[i, :]))
    end
    return W
end

function convolve(conv::Convolve, B, K)
    return _process_signal((X, _...) -> _convolve_impl(conv, X, K), B)
end

(conv::Convolve)(B, K) = convolve(conv, B, K)

function _lowpass_impl(conv::Convolve, B::AbstractMatrix{Float64}, scales::AbstractVector{Float64}, Bset, refactor::Bool, order::Int)
    _ = Bset
    W = Matrix{Float64}[]
    for scale in scales
        q = 1.0 / scale
        F = _factor!(conv, q; force=refactor || !haskey(conv.factor_cache, q))
        X = copy(B)
        for _ in 1:order
            X = F \ X
            X .*= q
        end
        push!(W, copy(X))
    end
    return W
end

function lowpass(conv::Convolve, B, scales=nothing; Bset=nothing, refactor::Bool=true, order::Int=1)
    scales === nothing && (scales = [1.0])
    return _process_signal((X, S, _...) -> _lowpass_impl(conv, X, S, Bset, refactor, order), B; scales=scales)
end

function _bandpass_impl(conv::Convolve, B::AbstractMatrix{Float64}, scales::AbstractVector{Float64}, order::Int)
    W = Matrix{Float64}[]
    for scale in scales
        q = 1.0 / scale
        F = _factor!(conv, q; force=true)
        inX = copy(B)
        for _ in 1:order
            X2 = F \ inX
            X1 = F \ X2
            inX = (4.0 * q) .* (conv.L * X1)
        end
        push!(W, copy(inX))
    end
    return W
end

function bandpass(conv::Convolve, B, scales=nothing; order::Int=1)
    scales === nothing && (scales = [1.0])
    return _process_signal((X, S, _...) -> _bandpass_impl(conv, X, S, order), B; scales=scales)
end

function _highpass_impl(conv::Convolve, B::AbstractMatrix{Float64}, scales::AbstractVector{Float64})
    W = Matrix{Float64}[]
    for scale in scales
        q = 1.0 / scale
        F = _factor!(conv, q; force=true)
        X1 = F \ B
        X2 = conv.L * X1
        push!(W, copy(X2))
    end
    return W
end

function highpass(conv::Convolve, B, scales=nothing)
    scales === nothing && (scales = [1.0])
    return _process_signal((X, S, _...) -> _highpass_impl(conv, X, S), B; scales=scales)
end

mutable struct DyConvolve
    n_vertices::Int
    L::SparseMatrixCSC{Float64, Int}
    poles::Vector{Float64}
    factors::Vector{SuiteSparse.CHOLMOD.Factor{Float64}}
    R::Any
    D::Vector{Float64}
end

function DyConvolve(L::SparseMatrixCSC, poles::AbstractVector)
    A = SparseMatrixCSC{Float64, Int}(L)
    p = Float64.(collect(poles))
    factors = [ldlt(Symmetric(A); shift=q) for q in p]
    return DyConvolve(size(A, 1), A, p, factors, nothing, Float64[])
end
DyConvolve(L::AbstractSparseMatrix, poles::AbstractVector) = DyConvolve(SparseMatrixCSC{Float64, Int}(L), poles)

function DyConvolve(L::SparseMatrixCSC, K::VFKernel)
    R, Q, D = _kernel_arrays(K)
    A = SparseMatrixCSC{Float64, Int}(L)
    factors = [ldlt(Symmetric(A); shift=q) for q in Q]
    return DyConvolve(size(A, 1), A, Float64.(Q), factors, R, D)
end
DyConvolve(L::AbstractSparseMatrix, K::VFKernel) = DyConvolve(SparseMatrixCSC{Float64, Int}(L), K)

function _dy_convolve_impl(conv::DyConvolve, B::AbstractMatrix{Float64})
    conv.R === nothing && throw(ArgumentError("Cannot call convolve without VFKernel poles/residues."))
    R = conv.R
    nDim = size(R, 2)
    W = zeros(Float64, size(B, 1), size(B, 2), nDim)
    if !isempty(conv.D)
        for d in eachindex(conv.D)
            W[:, :, d] .+= B .* conv.D[d]
        end
    end
    for (i, F) in enumerate(conv.factors)
        Z = F \ B
        _accum_residue!(W, Z, vec(R[i, :]))
    end
    return W
end

function convolve(conv::DyConvolve, B)
    return _process_signal((X, _...) -> _dy_convolve_impl(conv, X), B)
end

(conv::DyConvolve)(B) = convolve(conv, B)

function _dy_lowpass_impl(conv::DyConvolve, B::AbstractMatrix{Float64}, Bset, order::Int)
    _ = Bset
    W = Matrix{Float64}[]
    for (q, F) in zip(conv.poles, conv.factors)
        X = copy(B)
        for _ in 1:order
            X = F \ X
            X .*= q
        end
        push!(W, copy(X))
    end
    return W
end

function lowpass(conv::DyConvolve, B; Bset=nothing, order::Int=1)
    return _process_signal((X, _...) -> _dy_lowpass_impl(conv, X, Bset, order), B)
end

function _dy_bandpass_impl(conv::DyConvolve, B::AbstractMatrix{Float64}, order::Int)
    W = Matrix{Float64}[]
    for (q, F) in zip(conv.poles, conv.factors)
        inX = copy(B)
        for _ in 1:order
            X2 = F \ inX
            X1 = F \ X2
            inX = (4.0 * q) .* (conv.L * X1)
        end
        push!(W, copy(inX))
    end
    return W
end

function bandpass(conv::DyConvolve, B; order::Int=1)
    return _process_signal((X, _...) -> _dy_bandpass_impl(conv, X, order), B)
end

function _dy_highpass_impl(conv::DyConvolve, B::AbstractMatrix{Float64})
    W = Matrix{Float64}[]
    for F in conv.factors
        X1 = F \ B
        X2 = conv.L * X1
        push!(W, copy(X2))
    end
    return W
end

function highpass(conv::DyConvolve, B)
    return _process_signal((X, _...) -> _dy_highpass_impl(conv, X), B)
end

function addbranch(conv::DyConvolve, i::Int, j::Int, w::Real)
    if !(1 <= i <= conv.n_vertices && 1 <= j <= conv.n_vertices)
        return false
    end
    if w < 0
        throw(DomainError(w, "math domain error: weight w must be non-negative."))
    end

    ww = Float64(w)
    conv.L[i, j] -= ww
    conv.L[j, i] -= ww
    conv.L[i, i] += ww
    conv.L[j, j] += ww

    ws = sqrt(ww)
    C = sparsevec([i, j], [ws, -ws], conv.n_vertices)
    updown_sym = Symbol("updown!")
    if isdefined(SuiteSparse.CHOLMOD, updown_sym)
        updown! = getproperty(SuiteSparse.CHOLMOD, updown_sym)
        ok = true
        for F in conv.factors
            ok = ok && updown!(F, C, true)
        end
        return ok
    end

    # Fallback when rank-one updates are not exposed by this Julia build.
    conv.factors = [ldlt(Symmetric(conv.L); shift=q) for q in conv.poles]
    return true
end

addbranch!(conv::DyConvolve, i::Int, j::Int, w::Real) = addbranch(conv, i, j, w)
