using SpectralGraphWavelet
using Plots

include("filter_demo_plot_utils.jl")

function texas_bandpass_demo()
    L = resource("DELAY_TEXAS")
    C = resource("COORD_TEXAS")

    X = impulse(L, pyidx(600))
    s = 1e-1

    conv = Convolve(L)
    BP = bandpass(conv, X, s)
    p = plot_signal(BP, C; cmap="berlin", dot_size=3.0, dpi=400)
    return save_demo_plot(p, "texas_bandpass.png"; dpi=400)
end

function eastwest_bandpass_demo()
    L = resource("DELAY_EASTWEST")
    C = resource("COORD_EASTWEST")

    X = impulse(L, pyidx(65000))
    X .-= impulse(L, pyidx(35000))
    s = 3.0

    conv = Convolve(L)
    BP = bandpass(conv, X, s; order=3)
    p = plot_signal(BP, C; cmap="coolwarm", dot_size=3.0, dpi=400)
    return save_demo_plot(p, "eastwest_bandpass.png"; dpi=400)
end

function usa_multiscale_lowpass_demo()
    L = resource("DELAY_USA")
    C = resource("COORD_USA")

    X = 2.5 .* impulse(L, pyidx(22000))
    X .+= impulse(L, pyidx(72000))
    s = [0.5, 10.0, 100.0]

    conv = Convolve(L)
    Y = lowpass(conv, X, s; order=8)
    p = plot_signal(Y[1], C; cmap="viridis", dot_size=3.0, dpi=500)
    return save_demo_plot(p, "usa_multiscale_lowpass.png"; dpi=500)
end

function usa_bandpass_demo()
    L = resource("DELAY_USA")
    C = resource("COORD_USA")

    X = impulse(L, pyidx(75000))
    s = 2.0

    conv = Convolve(L)
    Y = bandpass(conv, X, s; order=10)
    p = plot_signal(Y, C; cmap="coolwarm", dot_size=3.0, dpi=500)
    return save_demo_plot(p, "usa_bandpass.png"; dpi=500)
end

function run_filter_demo_suite()
    outputs = String[]
    push!(outputs, texas_bandpass_demo())
    push!(outputs, eastwest_bandpass_demo())
    push!(outputs, usa_multiscale_lowpass_demo())
    push!(outputs, usa_bandpass_demo())
    return outputs
end

if abspath(PROGRAM_FILE) == @__FILE__
    outputs = run_filter_demo_suite()
    println("Generated filter demo plots:")
    for p in outputs
        println(" - ", p)
    end
end
