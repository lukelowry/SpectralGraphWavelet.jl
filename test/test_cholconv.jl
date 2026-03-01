@testset "cholconv" begin
    L = path_laplacian(80)
    x = impulse(L, 10)
    scales = [0.1, 1.0, 10.0]

    @testset "convolve analytical filters" begin
        conv = Convolve(L)
        lp = lowpass(conv, x, scales)
        bp = bandpass(conv, x, scales)
        hp = highpass(conv, x, scales)
        @test length(lp) == length(scales)
        @test length(bp) == length(scales)
        @test length(hp) == length(scales)
        @test all(length(v) == size(L, 1) for v in lp)
    end

    @testset "scalar scales and 1d behavior" begin
        conv = Convolve(L)
        y = lowpass(conv, x, 1.0)
        @test ndims(y) == 1
        z = highpass(conv, x, 1.0)
        @test ndims(z) == 1
    end

    @testset "complex signal handling" begin
        conv = Convolve(L)
        xc = complex.(x, x)
        yc = lowpass(conv, xc, 1.0)
        yr = lowpass(conv, x, 1.0)
        @test eltype(yc) <: Complex
        @test isapprox(yc, complex.(yr, yr); atol=1e-8)
    end

    @testset "vfkernel convolve" begin
        kernel_dict = resource("MODIFIED_MORLET")
        conv = Convolve(L)
        y = convolve(conv, x, kernel_dict)
        @test ndims(y) == 2
        @test size(y, 1) == size(L, 1)
        @test_throws ArgumentError convolve(conv, x, "not a kernel")

        bad = VFKernel(R=nothing, Q=nothing, D=nothing)
        @test_throws ArgumentError convolve(conv, x, bad)
    end

    @testset "direct term handling" begin
        conv = Convolve(L)
        K = VFKernel(R=[1.0;;], Q=[1.0], D=[5.0])
        y = convolve(conv, x, K)
        lp = lowpass(conv, x, 1.0)
        @test size(y, 2) == 1
        @test isapprox(vec(y[:, 1]), lp .+ 5.0 .* x; atol=1e-8)
    end

    @testset "dyconvolve filters and topology updates" begin
        poles = 1.0 ./ scales
        conv = DyConvolve(L, poles)
        lp_before = lowpass(conv, x)
        @test length(lp_before) == length(poles)
        ok = addbranch(conv, 10, 20, 1.0)
        @test ok
        lp_after = lowpass(conv, x)
        @test maximum(abs.(lp_before[1] .- lp_after[1])) > 1e-8
        @test !addbranch(conv, 0, 999, 1.0)
        @test_throws DomainError addbranch(conv, 10, 20, -1.0)
    end

    @testset "dyconvolve with VFKernel" begin
        vk = VFKernel(R=[1.0;;], Q=[1.0], D=[2.0])
        conv = DyConvolve(L, vk)
        y = convolve(conv, x)
        lp = lowpass(conv, x)[1]
        @test size(y, 2) == 1
        @test isapprox(vec(y[:, 1]), lp .+ 2.0 .* x; atol=1e-8)
    end
end
