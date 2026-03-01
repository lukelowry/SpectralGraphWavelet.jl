using Plots

const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "out")

"""
Convert a Python-style 0-based index to Julia's 1-based indexing.
"""
pyidx(i::Integer) = i + 1

function _colormap(name::AbstractString)
    lowered = lowercase(name)
    if lowered == "managua"
        # Tuned to match sgwt's Python "managua" visual appearance in docs:
        # diverging map with dark-plum center and warm/cool lobes.
        # Positive values map to cyan (as seen in usa_multiscale_lowpass_demo),
        # negative values map to orange (needed for usa_bandpass_demo).
        return cgrad(["#D89B59", "#5B2D57", "#4E5F9E", "#79CAE8"])
    end

    fallback = Dict(
        "berlin" => :RdBu,
        "coolwarm" => :coolwarm,
        "spectral" => :Spectral,
    )
    primary = Symbol(lowered)
    chosen = get(fallback, lowered, :viridis)

    try
        return cgrad(primary)
    catch
        return cgrad(chosen)
    end
end

function _xy_from_coords(coords::AbstractMatrix)
    c1 = vec(Float64.(coords[:, 1]))
    c2 = vec(Float64.(coords[:, 2]))

    c1_lon_like = minimum(c1) < -20.0
    c2_lon_like = minimum(c2) < -20.0

    if c1_lon_like && !c2_lon_like
        return c1, c2
    elseif c2_lon_like && !c1_lon_like
        return c2, c1
    end

    return c1, c2
end

function plot_signal(f::AbstractVecOrMat, coords::AbstractMatrix; cmap::AbstractString="Spectral", dot_size::Real=2.0, dpi::Int=400)
    fv = vec(Array(f))
    lon, lat = _xy_from_coords(coords)

    sabs = sort(abs.(fv))
    idx = max(1, length(sabs) - 19)
    mx = max(sabs[idx], eps(Float64))

    return scatter(
        lon,
        lat;
        marker_z=fv,
        color=_colormap(cmap),
        clims=(-mx, mx),
        markersize=dot_size,
        markerstrokewidth=0,
        legend=false,
        colorbar=false,
        showaxis=false,
        framestyle=:none,
        aspect_ratio=:equal,
        dpi=dpi,
    )
end

function save_demo_plot(p, filename::AbstractString; dpi::Int=400)
    mkpath(DEFAULT_OUTPUT_DIR)
    outpath = joinpath(DEFAULT_OUTPUT_DIR, filename)
    p.attr[:dpi] = dpi
    savefig(p, outpath)
    return outpath
end
