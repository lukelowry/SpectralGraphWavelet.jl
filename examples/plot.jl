"""
    plot_signal(f, S; cmap=nothing)

Parameters
    f: Signal data. Accepts Vector or (N x 1) Sparse/Dense Matrix.
    S: Coordinate Matrix (N x 2) [Longitude, Latitude]
    cmap: Colormap. Defaults to Red-Blue diverging gradient.
"""
function plot_signal(f::AbstractVecOrMat, S::AbstractMatrix; cmap=nothing)
    # 1. Convert to Dense 1D Vector
    f_vec = vec(Array(f))

    # 2. Extract Coordinates
    lon = S[:, 1]
    lat = S[:, 2]

    # 3. Robust Normalization
    sorted_abs = sort(abs.(f_vec))
    mx = sorted_abs[max(1, length(sorted_abs) - 19)]

    # 4. Handle Colormap
    if isnothing(cmap)
        # Wrap in cgrad() to create a Gradient from the palette
        color_scheme = cgrad(diverging_palette(0, 200, 32))
    else
        color_scheme = cmap
    end

    # 5. Create Plot
    p = scatter(lon, lat,
        zcolor = f_vec,
        c = color_scheme,
        clims = (-mx, mx),       
        markerstrokewidth = 0,   
        label = false,
        aspect_ratio = :equal,
        grid = false,
        title = "Signal Response",
        colorbar_ticks = nothing  # Removes the numeric labels from the colorbar
    )

    display(p)
    return p
end