"""
    lowpass(x; scale=1.0)

Analytical low-pass kernel response.
"""
function lowpass(x::AbstractArray; scale::Real=1.0)
    return @. 1.0 / (scale * x + 1.0)
end

"""
    highpass(x; scale=1.0)

Analytical high-pass kernel response.
"""
function highpass(x::AbstractArray; scale::Real=1.0)
    return @. (scale * x) / (scale * x + 1.0)
end

"""
    bandpass(x; scale=1.0, order=1)

Analytical band-pass wavelet kernel response.
"""
function bandpass(x::AbstractArray; scale::Real=1.0, order::Int=1)
    q = 1.0 / scale
    base = @. (4.0 * q * x) / (x + q)^2
    return base .^ order
end
