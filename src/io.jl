using MAT
using SparseArrays

"""
    load_laplacian(path::String)

Reads a .mat file containing the graph topology.
Expects key "A" (from your conversion script) or "L".
Returns: SparseMatrixCSC{Float64, Int}
"""
function load_laplacian(path::String; variable_name::String="A")
    if !isfile(path)
        error("File not found: $path")
    end

    vars = matread(path)

    # Check for the primary key "A", fallback to "L" if needed
    if !haskey(vars, variable_name)
        if haskey(vars, "L")
            variable_name = "L"
        else
            error("Matrix key '$variable_name' not found in $path. Keys found: $(keys(vars))")
        end
    end

    # Ensure correct type (Sparse Float64)
    return SparseMatrixCSC{Float64, Int}(vars[variable_name])
end

"""
    load_signal(path::String)

Reads a .mat file containing spatial signal data.
Expects keys "longitude" and "latitude".
Returns: (N x 2) Matrix{Float64} where column 1 is Longitude, column 2 is Latitude.
"""
function load_signal(path::String)
    if !isfile(path)
        error("File not found: $path")
    end

    data = matread(path)

    if !haskey(data, "longitude") || !haskey(data, "latitude")
        error("File $path must contain 'longitude' and 'latitude' keys.")
    end

    # vec() is crucial here: scipy.io.savemat often saves 1D arrays as 
    # (1 x N) or (N x 1) 2D matrices. vec() flattens them to simple vectors.
    lon = vec(data["longitude"])
    lat = vec(data["latitude"])

    if length(lon) != length(lat)
        error("Dimension mismatch in file $path: Longitude $(length(lon)) vs Latitude $(length(lat))")
    end

    # Combine into N x 2 Matrix
    return hcat(lon, lat)
end