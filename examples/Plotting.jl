using Plots
using SpectralGraphWavelet

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
        colorbar_ticks = nothing,  # Removes the numeric labels from the colorbar
        dpi = 400
    )

    savefig(raw"C:\Users\wyattluke.lowery\Documents\GitHub\SpectralGraphWavelet.jl\examples\USA_Wavelet_julia.png")
    return p
end

script_dir = dirname(@__FILE__)

path_L = joinpath(script_dir, "USA_DELAY.mat")
path_S = joinpath(script_dir, "usa_coords.mat")

A = load_laplacian(path_L)
S = load_signal(path_S)

scale = 10
conv = DyConvolve(A, [1/scale]) # pass poles (1/s)
b = impulse(A, 25000) 

c = bandpass(conv, b)[1]
c = bandpass(conv, c)[1]
c = bandpass(conv, c)[1]
c = bandpass(conv, c)[1]
results_bp = bandpass(conv, c)[1]

println("Attempting Plot")

f_vals = results_bp 
plot_signal(f_vals, S)

println("Plotted!")