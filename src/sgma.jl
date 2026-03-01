using Printf
using Statistics
import Base: close

const NetworkAnalysisResult = NamedTuple{(:peaks, :clusters), Tuple{Dict{String, Any}, Dict{String, Any}}}

mutable struct ModeTable
    frequency::Vector{Float64}
    damping::Vector{Float64}
    wavelength::Vector{Float64}
    magnitude::Vector{Float64}
    n_modes::Int
end

function ModeTable(frequency, damping, wavelength, magnitude)
    f = Float64.(vec(frequency))
    d = Float64.(vec(damping))
    w = Float64.(vec(wavelength))
    m = Float64.(vec(magnitude))
    return ModeTable(f, d, w, m, length(f))
end

function Base.show(io::IO, table::ModeTable)
    if table.n_modes == 0
        print(io, "ModeTable (empty - no modes identified)")
        return
    end
    println(io, "ModeTable ($(table.n_modes) mode$(table.n_modes > 1 ? "s" : "") identified)")
    println(io, repeat("-", 60))
    println(io, lpad("#", 3), "  ", lpad("Freq (Hz)", 10), "  ", lpad("Damping", 10), "  ", lpad("Wavelength", 12), "  ", lpad("Magnitude", 10))
    println(io, repeat("-", 60))
    for i in eachindex(table.frequency)
        println(
            io,
            lpad(string(i), 3),
            "  ",
            @sprintf("%10.4f", table.frequency[i]),
            "  ",
            @sprintf("%10.4f", table.damping[i]),
            "  ",
            @sprintf("%12.2f", table.wavelength[i]),
            "  ",
            @sprintf("%10.4f", table.magnitude[i]),
        )
    end
    print(io, repeat("-", 60))
end

to_dict(table::ModeTable) = Dict(
    "Frequency" => table.frequency,
    "Damping" => table.damping,
    "Wavelength" => table.wavelength,
    "Magnitude" => table.magnitude,
)

to_array(table::ModeTable) = hcat(table.frequency, table.damping, table.wavelength, table.magnitude)

mutable struct SGMA
    L::SparseMatrixCSC{Float64, Int}
    scales::Vector{Float64}
    freqs::Vector{Float64}
    order::Int
    w0::Float64
    Ts::Vector{Float64}
    wavlen::Vector{Float64}
    poles::Vector{Float64}
    _conv::Union{Nothing, DyConvolve}
    _B::Union{Nothing, Matrix{ComplexF64}}
    _t_cached::Union{Nothing, Vector{Float64}}
    _time_target_cached::Union{Nothing, Float64}
end

function SGMA(L::SparseMatrixCSC, scales::AbstractVector, freqs::AbstractVector; order::Int=10, w0::Real=2pi)
    A = SparseMatrixCSC{Float64, Int}(L)
    s = Float64.(collect(scales))
    f = Float64.(collect(freqs))
    Ts = w0 ./ (2pi .* f)
    wavlen = sqrt.(s)
    poles = 1.0 ./ s
    obj = SGMA(A, s, f, order, Float64(w0), Ts, wavlen, poles, nothing, nothing, nothing, nothing)
    finalizer(close, obj)
    return obj
end
SGMA(L::AbstractSparseMatrix, scales::AbstractVector, freqs::AbstractVector; order::Int=10, w0::Real=2pi) =
    SGMA(SparseMatrixCSC{Float64, Int}(L), scales, freqs; order=order, w0=w0)

function _get_conv(engine::SGMA)
    if engine._conv === nothing
        engine._conv = DyConvolve(engine.L, engine.poles)
    end
    return engine._conv
end

function _peak_local_max_fallback(image::AbstractMatrix{<:Real}; min_distance::Int=1, num_peaks::Int=typemax(Int))
    md = max(1, min_distance)
    nrow, ncol = size(image)
    local_coords = Tuple{Int, Int}[]
    for r in 1:nrow, c in 1:ncol
        image[r, c] > 0 || continue
        rlo = max(1, r - md)
        rhi = min(nrow, r + md)
        clo = max(1, c - md)
        chi = min(ncol, c + md)
        center = image[r, c]
        ismax = true
        for rr in rlo:rhi, cc in clo:chi
            if image[rr, cc] > center
                ismax = false
                break
            end
        end
        ismax && push!(local_coords, (r, c))
    end
    isempty(local_coords) && return zeros(Int, 0, 2)

    sort!(local_coords, by=rc -> image[rc[1], rc[2]], rev=true)
    suppressed = falses(nrow, ncol)
    final = Tuple{Int, Int}[]
    for (r, c) in local_coords
        suppressed[r, c] && continue
        push!(final, (r, c))
        length(final) == num_peaks && break
        rlo = max(1, r - md)
        rhi = min(nrow, r + md)
        clo = max(1, c - md)
        chi = min(ncol, c + md)
        suppressed[rlo:rhi, clo:chi] .= true
    end

    out = zeros(Int, length(final), 2)
    for (k, (r, c)) in enumerate(final)
        out[k, 1] = r
        out[k, 2] = c
    end
    return out
end

function _build_temporal_matrix(engine::SGMA, t::AbstractVector, time_target::Real)
    tf = Float64.(collect(t))
    if engine._B !== nothing &&
       engine._t_cached !== nothing &&
       length(tf) == length(engine._t_cached) &&
       all(isapprox.(tf, engine._t_cached; rtol=1e-8, atol=1e-8)) &&
       engine._time_target_cached == Float64(time_target)
        return engine._B
    end
    B = hcat([functions.gaussian_wavelet(tf; a=sc, b=time_target, w0=engine.w0) for sc in engine.Ts]...)
    engine._B = B
    engine._t_cached = copy(tf)
    engine._time_target_cached = Float64(time_target)
    return B
end

function spectrum(engine::SGMA, V::AbstractMatrix, t::AbstractVector; bus::Int, time::Real, VB=nothing, return_complex::Bool=false)
    n_nodes = size(engine.L, 1)
    if !(1 <= bus <= n_nodes)
        throw(ArgumentError("bus $bus out of bounds for $n_nodes nodes"))
    end

    VB_local = VB === nothing ? (V * _build_temporal_matrix(engine, t, time)) : VB
    conv = _get_conv(engine)
    spatial_responses = bandpass(conv, impulse(engine.L, bus); order=engine.order)
    A = hcat([vec(r) for r in spatial_responses]...)'
    Y = A * VB_local
    return return_complex ? Y : abs.(Y)
end

function analyze(engine::SGMA, V::AbstractMatrix, t::AbstractVector; bus::Int, time::Real, top_n::Int=5, min_dist::Int=5)
    return find_peaks(engine, spectrum(engine, V, t; bus=bus, time=time); top_n=top_n, min_dist=min_dist)
end

function find_peaks(engine::SGMA, spectrum::AbstractMatrix; top_n::Int=5, min_dist::Int=5, return_indices::Bool=false)
    Y = abs.(spectrum)
    coords = _peak_local_max_fallback(Y; min_distance=min_dist, num_peaks=top_n)
    keys = return_indices ? ["Wavelength", "Frequency", "Magnitude", "ScaleIdx", "FreqIdx"] : ["Wavelength", "Frequency", "Magnitude"]
    if size(coords, 1) == 0
        return Dict(k => Float64[] for k in keys)
    end

    mags = [Y[coords[i, 1], coords[i, 2]] for i in 1:size(coords, 1)]
    order = sortperm(mags; rev=true)
    sidx = [coords[i, 1] for i in order]
    fidx = [coords[i, 2] for i in order]
    result = Dict{String, Any}(
        "Wavelength" => engine.wavlen[sidx],
        "Frequency" => engine.freqs[fidx],
        "Magnitude" => mags[order],
    )
    if return_indices
        result["ScaleIdx"] = sidx
        result["FreqIdx"] = fidx
    end
    return result
end

function _gradient_axis(A::AbstractMatrix, x::AbstractVector, axis::Int)
    G = similar(A)
    n1, n2 = size(A)
    if axis == 1
        for j in 1:n2
            G[1, j] = (A[2, j] - A[1, j]) / (x[2] - x[1])
            for i in 2:(n1 - 1)
                G[i, j] = (A[i + 1, j] - A[i - 1, j]) / (x[i + 1] - x[i - 1])
            end
            G[n1, j] = (A[n1, j] - A[n1 - 1, j]) / (x[end] - x[end - 1])
        end
    else
        for i in 1:n1
            G[i, 1] = (A[i, 2] - A[i, 1]) / (x[2] - x[1])
            for j in 2:(n2 - 1)
                G[i, j] = (A[i, j + 1] - A[i, j - 1]) / (x[j + 1] - x[j - 1])
            end
            G[i, n2] = (A[i, n2] - A[i, n2 - 1]) / (x[end] - x[end - 1])
        end
    end
    return G
end

function find_modes(engine::SGMA, spectrum::AbstractMatrix; top_n::Int=5, min_dist::Int=5)
    if !(eltype(spectrum) <: Complex)
        throw(ArgumentError("find_modes requires complex spectrum (use return_complex=true)"))
    end
    Ymag = abs.(spectrum)
    peaks = find_peaks(engine, Ymag; top_n=top_n, min_dist=min_dist, return_indices=true)
    if isempty(peaks["Frequency"])
        return ModeTable(Float64[], Float64[], Float64[], Float64[])
    end

    log_s = log.(spectrum .+ 1e-20)
    grad_f = _gradient_axis(log_s, engine.freqs, 2)
    grad_s = _gradient_axis(log_s, engine.wavlen, 1)
    si = Int.(peaks["ScaleIdx"])
    fi = Int.(peaks["FreqIdx"])
    f0 = Float64.(peaks["Frequency"])
    s0 = Float64.(peaks["Wavelength"])
    omega_n = 2pi .* f0

    d_phi_df = [imag(grad_f[si[k], fi[k]]) for k in eachindex(si)]
    d_phi_ds = [imag(grad_s[si[k], fi[k]]) for k in eachindex(si)]

    zeta_f = @. -1.0 / (omega_n * d_phi_df)
    zeta_s = @. -s0 / (omega_n * d_phi_ds * s0^2 + 1e-10)
    wf = @. abs(d_phi_df) + 1e-10
    ws = @. abs(d_phi_ds * s0) + 1e-10
    damping = @. (wf * zeta_f + ws * zeta_s) / (wf + ws)

    singular = abs.(d_phi_df) .< 1e-8
    if any(singular)
        log_mag = log.(Ymag .+ 1e-20)
        d2 = _gradient_axis(_gradient_axis(log_mag, engine.freqs, 2), engine.freqs, 2)
        curv = [min(d2[si[k], fi[k]], -1e-10) for k in eachindex(si)]
        zeta_curv = @. sqrt(-2.0 / curv) / (2.0 * f0)
        for k in eachindex(damping)
            singular[k] && (damping[k] = zeta_curv[k])
        end
    end

    damping = [isfinite(z) ? clamp(z, 0.0, 1.0) : 0.0 for z in damping]
    return ModeTable(f0, damping, s0, peaks["Magnitude"])
end

function _compute_density_clusters(engine::SGMA, peaks::Dict{String, Any}, top_n::Int, min_dist::Int)
    if length(peaks["Wavelength"]) < 2
        return Dict("Wavelength" => Float64[], "Frequency" => Float64[], "Density" => Float64[])
    end
    try
        x = log10.(Float64.(peaks["Wavelength"]))
        y = Float64.(peaks["Frequency"])
        n = length(x)
        hx = max(1e-6, 1.06 * std(x) * n^(-1 / 5))
        hy = max(1e-6, 1.06 * std(y) * n^(-1 / 5))
        gx = log10.(engine.wavlen)
        gy = engine.freqs
        Z = zeros(Float64, length(gx), length(gy))
        c = 1.0 / (2pi * hx * hy * n)
        for i in eachindex(gx), j in eachindex(gy)
            acc = 0.0
            for k in 1:n
                dx = (gx[i] - x[k]) / hx
                dy = (gy[j] - y[k]) / hy
                acc += exp(-0.5 * (dx^2 + dy^2))
            end
            Z[i, j] = c * acc
        end
        peaks_z = find_peaks(engine, Z; top_n=top_n, min_dist=min_dist)
        return Dict(
            "Wavelength" => peaks_z["Wavelength"],
            "Frequency" => peaks_z["Frequency"],
            "Density" => peaks_z["Magnitude"],
        )
    catch
        return Dict("Wavelength" => Float64[], "Frequency" => Float64[], "Density" => Float64[])
    end
end

function analyze_many(engine::SGMA, V::AbstractMatrix, t::AbstractVector; time::Real, buses=nothing, top_n::Int=5, min_dist::Int=5, verbose::Bool=true)
    bus_list = buses === nothing ? collect(1:size(V, 1)) : collect(buses)
    VB = V * _build_temporal_matrix(engine, t, time)

    all_w = Vector{Float64}[]
    all_f = Vector{Float64}[]
    all_m = Vector{Float64}[]
    all_b = Vector{Int}[]
    for (idx, bus_idx) in enumerate(bus_list)
        p = find_peaks(engine, spectrum(engine, V, t; bus=bus_idx, time=time, VB=VB); top_n=top_n, min_dist=min_dist)
        if !isempty(p["Wavelength"])
            push!(all_w, Float64.(p["Wavelength"]))
            push!(all_f, Float64.(p["Frequency"]))
            push!(all_m, Float64.(p["Magnitude"]))
            push!(all_b, fill(bus_idx, length(p["Wavelength"])))
        end
        if verbose && idx % 50 == 0
            println("  Processed $idx/$(length(bus_list)) buses...")
        end
    end

    empty_peaks = Dict("Wavelength" => Float64[], "Frequency" => Float64[], "Magnitude" => Float64[], "Bus_ID" => Int[])
    empty_clusters = Dict("Wavelength" => Float64[], "Frequency" => Float64[], "Density" => Float64[])
    if isempty(all_w)
        return (peaks=empty_peaks, clusters=empty_clusters)
    end

    peaks = Dict{String, Any}(
        "Wavelength" => vcat(all_w...),
        "Frequency" => vcat(all_f...),
        "Magnitude" => vcat(all_m...),
        "Bus_ID" => vcat(all_b...),
    )
    clusters = _compute_density_clusters(engine, peaks, top_n, min_dist)
    return (peaks=peaks, clusters=clusters)
end

function close(engine::SGMA)
    engine._conv = nothing
    engine._B = nothing
    engine._t_cached = nothing
    engine._time_target_cached = nothing
    return nothing
end
