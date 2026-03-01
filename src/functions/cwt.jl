"""
    gaussian_wavelet(time; a=1.0, b=0.0, w0=1.0)

Morlet-like Gaussian wavelet used in SGMA.
"""
function gaussian_wavelet(time::AbstractVector{<:Real}; a::Real=1.0, b::Real=0.0, w0::Real=1.0)
    t = (Float64.(time) .- b) ./ a
    dt = length(time) > 1 ? (time[2] - time[1]) : 1.0
    norm_const = (dt / a) * pi^(-0.25)
    gauss = @. exp(-(t^2 / 2.0))
    ac = @. exp(im * w0 * t) - exp(-(w0^2) / 2.0)
    return norm_const .* gauss .* ac
end
