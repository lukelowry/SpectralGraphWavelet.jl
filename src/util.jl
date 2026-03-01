using MAT
using JSON
using LazyArtifacts
using Random
using SparseArrays
using Statistics

function _resolve_library_root()
    fallback = joinpath(@__DIR__, "..", "data", "library")
    isdir(fallback) && return fallback

    try
        root = artifact"sgwt_library"
        lib = joinpath(root, "library")
        return isdir(lib) ? lib : root
    catch
        throw(ArgumentError("Bundled library artifact 'sgwt_library' is unavailable."))
    end
end

const _LIBRARY_ROOT = _resolve_library_root()
const _RESOURCE_REGISTRY = Ref{Dict{String, Function}}(Dict{String, Function}())
const _RESOURCE_CACHE = Dict{String, Any}()

"""
    VFKernel

Rational vector-fitting kernel representation.
"""
struct VFKernel
    R::Any
    Q::Any
    D::Any
end

VFKernel(; R=nothing, Q=nothing, D=nothing) = VFKernel(R, Q, D)
VFKernel(data::AbstractDict) = VFKernel_from_dict(data)

# Backward-compatible alias with older Julia package versions.
const VFKern = VFKernel

"""
    ChebyKernel

Chebyshev polynomial approximation data for one or more kernels.
"""
struct ChebyKernel
    C::Matrix{Float64}
    spectrum_bound::Float64
    min_lambda::Float64
end

ChebyKernel(; C=zeros(0, 0), spectrum_bound=0.0, min_lambda=0.0) = ChebyKernel(
    Matrix{Float64}(C),
    Float64(spectrum_bound),
    Float64(min_lambda),
)
ChebyKernel(data::AbstractDict) = ChebyKernel_from_dict(data)

"""
    VFKernel_from_dict(data)

Construct a `VFKernel` from JSON-like dictionary data.
"""
function VFKernel_from_dict(data::AbstractDict)
    poles = get(data, "poles", Any[])
    rows = [vec(Float64.(get(p, "r", Any[]))) for p in poles]
    if isempty(rows)
        R = zeros(0, 0)
    else
        ncols = length(rows[1])
        any(length(r) != ncols for r in rows) &&
            throw(ArgumentError("All residue vectors must have the same length."))
        R = zeros(Float64, length(rows), ncols)
        for (i, r) in enumerate(rows)
            @inbounds R[i, :] = r
        end
    end
    Q = Float64[get(p, "q", 0.0) for p in poles]
    dval = get(data, "d", Float64[])
    D = dval isa Number ? Float64[dval] : Float64.(dval)
    return VFKernel(R=R, Q=Q, D=D)
end

"""
    ChebyKernel_from_dict(data)

Construct a `ChebyKernel` from dictionary data.
"""
function ChebyKernel_from_dict(data::AbstractDict)
    approximations = get(data, "approximations", Any[])
    bound = Float64(get(data, "spectrum_bound", 0.0))
    min_lambda = Float64(get(data, "min_lambda", 0.0))
    if isempty(approximations)
        return ChebyKernel(C=zeros(0, 0), spectrum_bound=bound, min_lambda=min_lambda)
    end

    coeffs = [Float64.(get(a, "coeffs", Any[])) for a in approximations]
    n = length(coeffs[1])
    if any(length(c) != n for c in coeffs)
        throw(ArgumentError("All 'coeffs' arrays must have the same length."))
    end
    C = hcat(coeffs...)
    return ChebyKernel(C=C, spectrum_bound=bound, min_lambda=min_lambda)
end

from_dict(::Type{VFKernel}, data::AbstractDict) = VFKernel_from_dict(data)
from_dict(::Type{ChebyKernel}, data::AbstractDict) = ChebyKernel_from_dict(data)

function _chebval_matrix(x::AbstractVector{<:Real}, C::AbstractMatrix{<:Real})
    n = size(C, 1)
    nd = size(C, 2)
    y = zeros(Float64, length(x), nd)
    if n == 0 || nd == 0
        return y
    end

    xv = Float64.(x)
    for d in 1:nd
        c = Float64.(C[:, d])
        b_kp1 = zeros(Float64, length(xv))
        b_kp2 = zeros(Float64, length(xv))
        for k in n:-1:2
            b_k = @. 2.0 * xv * b_kp1 - b_kp2 + c[k]
            b_kp2 = b_kp1
            b_kp1 = b_k
        end
        y[:, d] = @. xv * b_kp1 - b_kp2 + c[1]
    end
    return y
end

function _scale_x(kernel::ChebyKernel, x::AbstractVector{<:Real})
    return @. 2.0 * (x - kernel.min_lambda) / (kernel.spectrum_bound - kernel.min_lambda) - 1.0
end

"""
    evaluate(kernel::ChebyKernel, x)

Evaluate a Chebyshev kernel at points `x`.
"""
function evaluate(kernel::ChebyKernel, x::AbstractVector{<:Real})
    if isempty(kernel.C)
        return zeros(length(x), 0)
    end
    y = _chebval_matrix(_scale_x(kernel, x), kernel.C)
    return size(y, 2) == 1 ? vec(y) : y
end

"""
    impulse(L, n=1, n_timesteps=1)

Generate a Dirac impulse on vertex `n`.
"""
function impulse(L::AbstractSparseMatrix, n::Int=1, n_timesteps::Int=1)
    if n_timesteps == 1
        b = zeros(Float64, size(L, 1))
        b[n] = 1.0
        return b
    end

    b = zeros(Float64, size(L, 1), n_timesteps)
    b[n, :] .= 1.0
    return b
end

"""
    estimate_spectral_bound(L)

Estimate the largest eigenvalue using power iteration and a 1.01 margin.
"""
function estimate_spectral_bound(L::AbstractMatrix; maxiter::Int=200, tol::Float64=1e-8)
    n = size(L, 1)
    x = randn(n)
    x ./= norm(x)
    lambda_prev = 0.0
    lambda_est = 0.0
    for _ in 1:maxiter
        y = L * x
        ny = norm(y)
        ny == 0 && return 0.0
        x = y ./ ny
        lambda_est = dot(x, L * x)
        if abs(lambda_est - lambda_prev) <= tol * max(1.0, abs(lambda_est))
            break
        end
        lambda_prev = lambda_est
    end
    return Float64(abs(lambda_est) * 1.01)
end

function _load_resource(path::AbstractString, loader::Function)
    if !isfile(path)
        throw(ArgumentError("Resource not found: $path"))
    end
    return loader(path)
end

function _mat_loader(path::AbstractString; to_csc::Bool=false)
    data = matread(path)
    vars = [k for k in keys(data) if !startswith(String(k), "__")]
    isempty(vars) && throw(ArgumentError("No data variables found in MAT file: $path"))

    raw = data[first(vars)]
    if to_csc
        if raw isa SparseMatrixCSC
            return SparseMatrixCSC{Float64, Int}(raw)
        end
        return sparse(Float64.(Array(raw)))
    end

    if length(vars) > 1
        cols = [vec(Float64.(Array(data[k]))) for k in vars]
        return hcat(cols...)
    end

    arr = Array(raw)
    if ndims(arr) == 2 && size(arr, 1) == 1
        return permutedims(arr)
    end
    return arr
end

function _json_kernel_loader(path::AbstractString)
    open(path, "r") do io
        return JSON.parse(read(io, String))
    end
end

const _PLY_TYPE_READERS = Dict{String, DataType}(
    "char" => Int8,
    "uchar" => UInt8,
    "short" => Int16,
    "ushort" => UInt16,
    "int" => Int32,
    "uint" => UInt32,
    "float" => Float32,
    "double" => Float64,
)

function _read_ply_scalar(io::IO, typename::String)
    T = get(_PLY_TYPE_READERS, typename, Float32)
    return read(io, T)
end

function _parse_ply(filepath::AbstractString)
    open(filepath, "r") do io
        fmt = "ascii"
        vertex_count = 0
        face_count = 0
        vertex_props = Tuple{String, String}[]
        current_element = ""

        while !eof(io)
            line = strip(readline(io))
            line == "end_header" && break
            isempty(line) && continue
            parts = split(line)
            isempty(parts) && continue
            if parts[1] == "format"
                fmt = parts[2]
            elseif parts[1] == "element"
                current_element = parts[2]
                if current_element == "vertex"
                    vertex_count = parse(Int, parts[3])
                elseif current_element == "face"
                    face_count = parse(Int, parts[3])
                end
            elseif parts[1] == "property" && current_element == "vertex"
                # property <type> <name>
                push!(vertex_props, (parts[3], parts[2]))
            end
        end

        vertices = NTuple{3, Float64}[]
        faces = Vector{Int}[]

        if fmt == "ascii"
            for _ in 1:vertex_count
                vals = split(strip(readline(io)))
                push!(vertices, (parse(Float64, vals[1]), parse(Float64, vals[2]), parse(Float64, vals[3])))
            end
            for _ in 1:face_count
                vals = split(strip(readline(io)))
                n = parse(Int, vals[1])
                idx = [parse(Int, vals[k + 1]) for k in 1:n]
                push!(faces, idx)
            end
        elseif fmt == "binary_little_endian"
            prop_names = [p[1] for p in vertex_props]
            x_idx = findfirst(==("x"), prop_names)
            y_idx = findfirst(==("y"), prop_names)
            z_idx = findfirst(==("z"), prop_names)
            xyz_idx = if x_idx !== nothing && y_idx !== nothing && z_idx !== nothing
                (x_idx, y_idx, z_idx)
            else
                (1, 2, 3)
            end

            for _ in 1:vertex_count
                vals = Float64[]
                for (_, ptype) in vertex_props
                    push!(vals, Float64(_read_ply_scalar(io, ptype)))
                end
                push!(vertices, (vals[xyz_idx[1]], vals[xyz_idx[2]], vals[xyz_idx[3]]))
            end

            for _ in 1:face_count
                n = Int(read(io, UInt8))
                idx = Int[]
                for _ in 1:n
                    push!(idx, Int(read(io, Int32)))
                end
                push!(faces, idx)
            end
        else
            throw(ArgumentError("Unsupported PLY format: $fmt"))
        end

        return vertices, faces, vertex_count
    end
end

"""
    load_ply_laplacian(filepath)

Load a PLY mesh and return its graph Laplacian `L = B*B'`.
"""
function load_ply_laplacian(filepath::AbstractString)
    vertices, faces, vertex_count = _parse_ply(filepath)
    _ = vertices

    unique_edges = Set{Tuple{Int, Int}}()
    for face in faces
        n = length(face)
        for i in 1:n
            u0 = face[i]
            v0 = face[mod1(i + 1, n)]
            u = min(u0, v0)
            v = max(u0, v0)
            push!(unique_edges, (u + 1, v + 1)) # convert 0-based -> 1-based
        end
    end

    edges = collect(unique_edges)
    sort!(edges)
    num_edges = length(edges)

    rows = Int[]
    cols = Int[]
    vals = Float64[]
    sizehint!(rows, 2 * num_edges)
    sizehint!(cols, 2 * num_edges)
    sizehint!(vals, 2 * num_edges)

    for (k, (u, v)) in enumerate(edges)
        push!(rows, u); push!(cols, k); push!(vals, 1.0)
        push!(rows, v); push!(cols, k); push!(vals, -1.0)
    end

    B = sparse(rows, cols, vals, vertex_count, num_edges)
    return B * B'
end

"""
    load_ply_xyz(filepath)

Load vertex coordinates from a PLY file as `(N, 3)` matrix.
"""
function load_ply_xyz(filepath::AbstractString)
    vertices, _, _ = _parse_ply(filepath)
    xyz = zeros(Float64, length(vertices), 3)
    for (i, (x, y, z)) in enumerate(vertices)
        xyz[i, 1] = x
        xyz[i, 2] = y
        xyz[i, 3] = z
    end
    return xyz
end

function load_laplacian(path::AbstractString; variable_name::String="A")
    if !isfile(path)
        throw(ArgumentError("File not found: $path"))
    end
    vars = matread(path)
    key = variable_name
    if !haskey(vars, key)
        key = haskey(vars, "L") ? "L" : ""
    end
    isempty(key) && throw(ArgumentError("Matrix key '$variable_name' not found in $path. Keys: $(collect(keys(vars)))"))
    return sparse(Float64.(Array(vars[key])))
end

function load_signal(path::AbstractString)
    if !isfile(path)
        throw(ArgumentError("File not found: $path"))
    end
    data = matread(path)
    if !haskey(data, "longitude") || !haskey(data, "latitude")
        throw(ArgumentError("File $path must contain 'longitude' and 'latitude' keys."))
    end
    lon = vec(Float64.(Array(data["longitude"])))
    lat = vec(Float64.(Array(data["latitude"])))
    length(lon) == length(lat) || throw(ArgumentError("Dimension mismatch in $path"))
    return hcat(lon, lat)
end

function _lap(name::String, region::String)
    path = joinpath(_LIBRARY_ROOT, name, "$region.mat")
    return _load_resource(path, p -> _mat_loader(p; to_csc=true))
end

function _sig(region::String)
    path = joinpath(_LIBRARY_ROOT, "SIGNALS", "$region.mat")
    return _load_resource(path, p -> _mat_loader(p; to_csc=false))
end

function _kern(name::String)
    path = joinpath(_LIBRARY_ROOT, "KERNELS", "$name.json")
    return _load_resource(path, _json_kernel_loader)
end

function _ply_lap(name::String)
    path = joinpath(_LIBRARY_ROOT, "MESH", "$name.ply")
    return _load_resource(path, load_ply_laplacian)
end

function _ply_xyz(name::String)
    path = joinpath(_LIBRARY_ROOT, "MESH", "$name.ply")
    return _load_resource(path, load_ply_xyz)
end

function _discover_resources()
    registry = Dict{String, Function}()
    if !isdir(_LIBRARY_ROOT)
        return registry
    end

    folder_configs = [
        ("KERNELS", ".json", s -> s, s -> (() -> _kern(s)), nothing),
        ("DELAY", ".mat", s -> "DELAY_$s", s -> (() -> _lap("DELAY", s)), nothing),
        ("IMPEDANCE", ".mat", s -> "IMPEDANCE_$s", s -> (() -> _lap("IMPEDANCE", s)), nothing),
        ("LENGTH", ".mat", s -> "LENGTH_$s", s -> (() -> _lap("LENGTH", s)), nothing),
        ("SIGNALS", ".mat", s -> "COORD_$s", s -> (() -> _sig(s)), nothing),
        (
            "MESH",
            ".ply",
            s -> "MESH_$s",
            s -> (() -> _ply_lap(s)),
            s -> ("$(s)_XYZ", () -> _ply_xyz(s)),
        ),
    ]

    for (folder, ext, key_fn, loader_fn, secondary_fn) in folder_configs
        dirpath = joinpath(_LIBRARY_ROOT, folder)
        isdir(dirpath) || continue
        for item in readdir(dirpath)
            endswith(item, ext) || continue
            stem = item[1:(end - length(ext))]
            registry[key_fn(stem)] = loader_fn(stem)
            if secondary_fn !== nothing
                key, loader = secondary_fn(stem)
                registry[key] = loader
            end
        end
    end

    return registry
end

function _ensure_registry()
    if isempty(_RESOURCE_REGISTRY[])
        _RESOURCE_REGISTRY[] = _discover_resources()
    end
    return _RESOURCE_REGISTRY[]
end

function _get_resource(name::String)
    if haskey(_RESOURCE_CACHE, name)
        return _RESOURCE_CACHE[name]
    end
    registry = _ensure_registry()
    haskey(registry, name) || throw(ArgumentError("Unknown resource: $name"))
    value = registry[name]()
    _RESOURCE_CACHE[name] = value
    return value
end

function _resource_names()
    return sort!(collect(keys(_ensure_registry())))
end

"""
    resource(name)

Load a bundled resource by name.
"""
resource(name::AbstractString) = _get_resource(String(name))
resource(name::Symbol) = _get_resource(String(name))

"""
    resource_names()

Return sorted bundled resource names.
"""
resource_names() = _resource_names()

"""
    list_graphs()

Print available bundled graph resources.
"""
function list_graphs(io::IO=stdout)
    registry = _ensure_registry()
    graph_prefixes = ("DELAY_", "IMPEDANCE_", "LENGTH_", "MESH_")
    graph_keys = [k for k in keys(registry) if any(startswith(k, p) for p in graph_prefixes)]
    if isempty(graph_keys)
        println(io, "No graphs found in the library.")
        return nothing
    end

    println(io, rpad("Graph Name", 30), rpad("Vertices", 10), "Edges")
    println(io, repeat("-", 52))
    for key in sort(graph_keys)
        try
            L = _get_resource(key)
            if L isa AbstractMatrix && hasproperty(L, :nnz)
                n_vertices = size(L, 1)
                diag_nnz = count(!iszero, diag(L))
                n_edges = Int(div(nnz(L) - diag_nnz, 2))
                println(io, rpad(key, 30), rpad(string(n_vertices), 10), n_edges)
            end
        catch
            continue
        end
    end
    return nothing
end

"""
    get_cholmod_dll()

Compatibility stub for Python API parity.
"""
function get_cholmod_dll()
    throw(ArgumentError("get_cholmod_dll is not needed in Julia; native SuiteSparse is used."))
end

"""
    get_klu_dll()

Compatibility stub for Python API parity.
"""
function get_klu_dll()
    throw(ArgumentError("get_klu_dll is not needed in Julia; native SuiteSparse is used."))
end
