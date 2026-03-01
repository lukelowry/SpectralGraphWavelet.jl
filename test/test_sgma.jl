@testset "sgma" begin
    L = path_laplacian(20)
    scales = exp10.(range(-1, 1, length=5))
    freqs = collect(range(0.1, 1.0, length=5))
    engine = SGMA(L, scales, freqs)
    V = random_signal(L; ncols=12)
    t = collect(range(0.0, 5.0, length=size(V, 2)))

    @testset "initialization" begin
        @test length(engine.scales) == 5
        @test length(engine.freqs) == 5
        @test length(engine.Ts) == 5
        @test engine._conv === nothing
    end

    @testset "spectrum shape and caching" begin
        Y = spectrum(engine, V, t; bus=1, time=2.5)
        @test size(Y) == (length(scales), length(freqs))
        @test all(Y .>= 0)

        B1 = SpectralGraphWavelet._build_temporal_matrix(engine, t, 2.5)
        B2 = SpectralGraphWavelet._build_temporal_matrix(engine, t, 2.5)
        @test B1 === B2
    end

    @testset "find_peaks and analyze" begin
        Y = zeros(5, 5)
        Y[3, 3] = 10.0
        p = find_peaks(engine, Y; top_n=1, min_dist=1)
        @test !isempty(p["Wavelength"])
        @test p["Magnitude"][1] == 10.0

        pa = analyze(engine, V, t; bus=1, time=2.5, top_n=2, min_dist=1)
        @test haskey(pa, "Frequency")
    end

    @testset "analyze_many result shape" begin
        res = analyze_many(engine, V, t; time=2.5, buses=collect(1:10), verbose=false, min_dist=1)
        @test haskey(res.peaks, "Bus_ID")
        @test haskey(res.clusters, "Density")
    end

    @testset "invalid bus" begin
        @test_throws ArgumentError spectrum(engine, V, t; bus=100, time=2.5)
    end

    @testset "mode table and find_modes" begin
        Yc = spectrum(engine, V, t; bus=1, time=2.5, return_complex=true)
        modes = find_modes(engine, Yc; top_n=3, min_dist=1)
        @test modes isa ModeTable
        @test modes.n_modes >= 0
        d = to_dict(modes)
        @test haskey(d, "Frequency")
        arr = to_array(modes)
        @test size(arr, 2) == 4

        @test_throws ArgumentError find_modes(engine, abs.(Yc))
    end

    close(engine)
end
