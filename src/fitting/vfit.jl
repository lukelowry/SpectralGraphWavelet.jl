using LinearAlgebra
using Statistics

function _cauchy(s::AbstractVector{<:Complex}, poles::AbstractVector{<:Complex})
    return 1.0 ./ (s .- permutedims(poles))
end

function _lstsq(A::AbstractMatrix, b::AbstractArray)
    return A \ b
end

function _basis(s::AbstractVector{<:Complex}, poles::AbstractVector{<:Complex}, fit_d::Bool, fit_e::Bool)
    M = Matrix{ComplexF64}(_cauchy(s, poles))
    if fit_d
        M = hcat(M, ones(ComplexF64, length(s), 1))
    end
    if fit_e
        M = hcat(M, reshape(ComplexF64.(s), :, 1))
    end
    return M
end

function _enforce_conjugates(poles::AbstractVector{<:Complex})
    N = length(poles)
    tol = 1e-6 * (maximum(abs.(poles)) + 1.0)
    realp = ComplexF64[real(z) + 0im for z in poles if abs(imag(z)) < tol]
    upper = ComplexF64[z for z in poles if imag(z) >= tol]
    if isempty(upper) && !isempty(poles)
        upper = ComplexF64[conj(z) for z in poles if imag(z) <= -tol]
    end
    allp = vcat(realp, upper, conj.(upper))
    if isempty(allp)
        allp = ComplexF64.(poles)
    end
    sort!(allp, by=x -> (real(x), imag(x)))
    if length(allp) > N
        allp = allp[1:N]
    elseif length(allp) < N
        fallback = sort(ComplexF64.(poles), by=x -> (real(x), imag(x)))
        append!(allp, fallback[1:(N - length(allp))])
    end
    return allp
end

function _initial_poles(s::AbstractVector{<:Complex}, N::Int, real_only::Bool)
    _geomspace(lo::Float64, hi::Float64, n::Int) = n <= 1 ? [lo] : exp10.(range(log10(lo), log10(hi), length=n))

    if real_only || N < 2
        r = abs.(s[s .!= 0])
        if isempty(r)
            r = abs.(s) .+ 1
        end
        lo = max(minimum(r) * 0.1, 1e-6)
        hi = maximum(r) * 2.0
        vals = _geomspace(lo, hi, N)
        return ComplexF64.(-vals)
    end

    w = abs.(imag.(s))
    w = w[w .> 0]
    if isempty(w)
        return ComplexF64.(-_geomspace(1e-2, 1e2, N))
    end
    lo = max(minimum(w) * 0.5, 1e-3)
    hi = maximum(w) * 1.5
    pts = _geomspace(lo, hi, max(1, div(N, 2)))
    cc = ComplexF64[Complex(-0.01 * p, p) for p in pts]
    poles = vcat(cc, conj.(cc))
    if isodd(N)
        poles = vcat(poles, ComplexF64(-pts[end], 0.0))
    end
    return poles[1:N]
end

function _relocate(
    s::AbstractVector{<:Complex},
    f::AbstractMatrix{<:Complex},
    poles::AbstractVector{<:Complex},
    fit_d::Bool,
    fit_e::Bool,
    real_only::Bool,
)
    K, A = size(f)
    N = length(poles)

    Phi_sig = _basis(s, poles, true, false)
    Phi = _basis(s, poles, fit_d, fit_e)
    qrf = qr(Phi)
    Q = Matrix(qrf.Q)[:, 1:size(Phi, 2)]

    weighted = Array{ComplexF64}(undef, K, A, N + 1)
    for i in 1:K, a in 1:A, j in 1:(N + 1)
        weighted[i, a, j] = f[i, a] * Phi_sig[i, j]
    end

    # Match NumPy's row-major reshape(K, -1): columns index (a, j) with j fastest.
    flat = Matrix{ComplexF64}(undef, K, A * (N + 1))
    col = 1
    for a in 1:A
        for j in 1:(N + 1)
            @views flat[:, col] .= weighted[:, a, j]
            col += 1
        end
    end

    proj = Q * (Q' * flat)
    rows = Matrix{ComplexF64}(undef, K * A, N + 1)
    for i in 1:K
        for a in 1:A
            r = (i - 1) * A + a
            base = (a - 1) * (N + 1)
            for j in 1:(N + 1)
                rows[r, j] = flat[i, base + j] - proj[i, base + j]
            end
        end
    end

    w = sqrt(K * A)
    M = vcat(rows, reshape(sum(Phi_sig; dims=1) .* w, 1, :))
    rhs = vcat(zeros(ComplexF64, size(M, 1) - 1), ComplexF64(K * w))
    z = _lstsq(M, rhs)
    c_tilde = z[1:N]
    d_tilde = z[N + 1] + 1e-30

    H = Diagonal(poles) - ones(ComplexF64, N, 1) * reshape(c_tilde ./ d_tilde, 1, :)
    new_poles = eigvals(H)
    new_poles = ComplexF64.(-abs.(real.(new_poles)), imag.(new_poles))

    if real_only
        rp = ComplexF64.(real.(new_poles))
        sort!(rp, by=real)
        return rp
    end
    return _enforce_conjugates(new_poles)
end

mutable struct VFModel
    poles::Vector{ComplexF64}
    residues::Matrix{ComplexF64}
    d::Vector{ComplexF64}
    e::Vector{ComplexF64}
    rmse::Float64
    iters::Int
end

function Base.show(io::IO, m::VFModel)
    print(
        io,
        "VFModel(",
        length(m.poles),
        " poles, ",
        size(m.residues, 2),
        " ch, rmse=",
        @sprintf("%.3e", m.rmse),
        ", ",
        m.iters,
        " iters)",
    )
end

function (m::VFModel)(s)
    sv = ComplexF64.(vec(s))
    return _cauchy(sv, m.poles) * m.residues .+ permutedims(m.d) .+ (reshape(sv, :, 1) .* permutedims(m.e))
end

function fit(
    ::Type{VFModel},
    s,
    f,
    n_poles::Int;
    poles=nothing,
    iters::Int=30,
    tol::Float64=1e-9,
    fit_d::Bool=true,
    fit_e::Bool=false,
    real_only::Bool=false,
)
    sv = ComplexF64.(vec(s))
    fv = ComplexF64.(f)
    F = ndims(fv) == 1 ? reshape(fv, :, 1) : fv
    if size(F, 1) != length(sv)
        F = permutedims(F)
    end

    a = poles === nothing ? _initial_poles(sv, n_poles, real_only) : ComplexF64.(poles)
    it_used = 0
    for it in 1:iters
        a_prev = copy(a)
        a = _relocate(sv, F, a, fit_d, fit_e, real_only)
        it_used = it
        rel = maximum(abs.(a .- a_prev) ./ (abs.(a_prev) .+ 1e-15))
        rel < tol && break
    end

    X = _lstsq(_basis(sv, a, fit_d, fit_e), F)
    N = length(a)
    A = size(F, 2)
    d = fit_d ? vec(X[N + 1, :]) : zeros(ComplexF64, A)
    e = fit_e ? vec(X[N + 1 + (fit_d ? 1 : 0), :]) : zeros(ComplexF64, A)
    residues = X[1:N, :]

    model = VFModel(a, residues, d, e, 0.0, it_used)
    err = abs.(F .- model(sv))
    model.rmse = sqrt(mean(err .^ 2))
    return model
end

function kernel(::Type{VFModel}, lam, g, n_poles::Int; kwargs...)
    kw = Dict{Symbol, Any}(kwargs)
    kw[:fit_d] = get(kw, :fit_d, true)
    kw[:fit_e] = false
    model = fit(VFModel, Float64.(lam), Float64.(g), n_poles; real_only=true, kw...)
    return VFKernel(R=real.(model.residues), Q=real.(model.poles), D=real.(model.d))
end

function Base.getproperty(::Type{VFModel}, s::Symbol)
    if s === :fit
        return (args...; kwargs...) -> fit(VFModel, args...; kwargs...)
    elseif s === :kernel
        return (args...; kwargs...) -> kernel(VFModel, args...; kwargs...)
    end
    return getfield(VFModel, s)
end
