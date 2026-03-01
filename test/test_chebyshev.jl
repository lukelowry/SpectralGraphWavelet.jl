@testset "chebyconv" begin
    L = path_laplacian(30)
    I_sig = Matrix(I, size(L, 1), size(L, 1))
    R_sig = random_signal(L; ncols=4)

    @testset "single coefficient kernel" begin
        ub = estimate_spectral_bound(L)
        K = ChebyKernel(C=[3.0;;], spectrum_bound=ub, min_lambda=0.0)
        conv = ChebyConvolve(L)
        y = convolve(conv, I_sig, K)
        @test isapprox(y[:, :, 1], 3.0 .* I_sig; atol=1e-8)
    end

    @testset "identity kernel returns input" begin
        ub = estimate_spectral_bound(L)
        C = zeros(2, 1)
        C[1, 1] = 1.0
        K = ChebyKernel(C=C, spectrum_bound=ub, min_lambda=0.0)
        conv = ChebyConvolve(L)
        y = convolve(conv, I_sig, K)
        @test isapprox(y[:, :, 1], I_sig; atol=1e-8)
    end

    @testset "high-order stability" begin
        f = x -> exp.(-x)
        K = ChebyModel.kernel(L, f, 30)
        conv = ChebyConvolve(L)
        y = convolve(conv, I_sig, K)
        @test !any(isnan, y)
        @test !any(isinf, y)
    end

    @testset "1d and complex input" begin
        f = x -> exp.(-x)
        K = ChebyModel.kernel(L, f, 10)
        conv = ChebyConvolve(L)
        v = ones(size(L, 1))
        y = convolve(conv, v, K)
        @test ndims(y) == 2

        vc = complex.(v, v)
        yc = convolve(conv, vc, K)
        @test eltype(yc) <: Complex
    end

    @testset "convolve_multi" begin
        K1 = ChebyModel.kernel(L, (x -> exp.(-x)), 8)
        K2 = ChebyModel.kernel(L, (x -> @. 1.0 / (x + 1.0)), 8)
        conv = ChebyConvolve(L)
        ys = convolve_multi(conv, R_sig, [K1, K2])
        @test length(ys) == 2
        @test isapprox(ys[1], convolve(conv, R_sig, K1); atol=1e-8)
        @test isapprox(ys[2], convolve(conv, R_sig, K2); atol=1e-8)
    end
end
