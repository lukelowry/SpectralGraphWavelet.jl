@testset "functions" begin
    x = collect(range(0.0, 10.0, length=200))

    @testset "lowpass" begin
        y = SpectralGraphWavelet.functions.lowpass(x; scale=1.0)
        @test isapprox(y[1], 1.0)
        @test all(diff(y) .<= 1e-12)
    end

    @testset "highpass" begin
        y = SpectralGraphWavelet.functions.highpass(x; scale=1.0)
        @test isapprox(y[1], 0.0)
        @test all(diff(y) .>= -1e-12)
    end

    @testset "bandpass" begin
        y = SpectralGraphWavelet.functions.bandpass(x; scale=1.0)
        @test isapprox(y[1], 0.0)
        peak = x[argmax(y)]
        @test abs(peak - 1.0) < 0.2
    end

    @testset "gaussian wavelet" begin
        t = collect(range(0.0, 10.0, length=1000))
        psi = SpectralGraphWavelet.functions.gaussian_wavelet(t; a=0.5, b=5.0, w0=2pi)
        @test length(psi) == length(t)
        @test eltype(psi) <: Complex
    end
end
