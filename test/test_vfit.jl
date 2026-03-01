@testset "vfit" begin
    @testset "fit recovers rational response" begin
        s = im .* collect(range(1.0, 100.0, length=300))
        f = @. 5.0 + 0.1 * s + 1.0 / (s + 3.0) + 2.0 / (s + 20.0)
        m = VFModel.fit(s, f, 2; fit_d=true, fit_e=true)
        rp = sort(real.(m.poles))
        @test abs(rp[1] + 20.0) < 2.0
        @test abs(rp[2] + 3.0) < 2.0
        @test abs(real(m.d[1]) - 5.0) < 1.0
        @test abs(real(m.e[1]) - 0.1) < 0.05
        @test m.rmse < 1e-2
    end

    @testset "fit multi-channel" begin
        s = im .* collect(range(1.0, 50.0, length=200))
        f = hcat(@.(1.0 / (s + 3.0)), @.(2.0 / (s + 10.0)))
        m = VFModel.fit(s, f, 2)
        @test size(m.residues, 2) == 2
        @test length(m.d) == 2
        @test m.rmse < 1e-1
    end

    @testset "kernel returns constrained VFKernel" begin
        lam = collect(range(0.01, 100.0, length=300))
        g = @. 1.0 / (lam + 1.0)
        K = VFModel.kernel(lam, g, 4)
        @test K isa VFKernel
        @test all(abs.(imag.(K.Q)) .< 1e-10)
        @test all(real.(K.Q) .< 0.0)
        @test size(K.R, 1) == 4
    end
end
