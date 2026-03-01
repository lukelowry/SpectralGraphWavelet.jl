@testset "chebymodel" begin
    L = path_laplacian(40)

    @testset "fit approximates function" begin
        f = x -> @. 1.0 / (x + 1.0)
        m = ChebyModel.fit(f, 40, 100.0)
        x = collect(range(0.01, 100.0, length=200))
        err = maximum(abs.(m(x) .- f(x)))
        @test err < 5e-2
    end

    @testset "kernel returns ChebyKernel" begin
        f = x -> @. 1.0 / (x + 1.0)
        K = ChebyModel.kernel(L, f, 20)
        @test K isa ChebyKernel
        @test K.spectrum_bound > 0.01
    end

    @testset "invalid order" begin
        f = x -> x
        @test_throws ArgumentError ChebyModel.kernel(L, f, 0)
    end

    @testset "multi-output function" begin
        f = x -> hcat(exp.(-x), sin.(x))
        K = ChebyModel.kernel(L, f, 8)
        @test size(K.C, 2) == 2
        y = evaluate(K, collect(range(0.0, K.spectrum_bound, length=10)))
        @test size(y) == (10, 2)
    end

    @testset "sampling strategies" begin
        f = x -> exp.(-x)
        for s in ("linear", "quadratic", "logarithmic")
            K = ChebyModel.kernel(L, f, 6; sampling=s)
            @test size(K.C, 1) >= 1
        end
    end

    @testset "adaptive fit" begin
        f = x -> exp.(-x)
        K = ChebyModel.kernel(L, f, 6; adaptive=true, target_error=0.01, max_order=40)
        @test size(K.C, 1) >= 1
    end
end
