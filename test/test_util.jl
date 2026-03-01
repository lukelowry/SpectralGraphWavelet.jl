@testset "util" begin
    @testset "vfkernel parsing" begin
        d = Dict(
            "poles" => Any[
                Dict("q" => 1.0, "r" => [0.1, 0.2]),
                Dict("q" => 2.0, "r" => [0.3, 0.4]),
            ],
            "d" => [0.5, 0.6],
        )
        k = from_dict(VFKernel, d)
        @test k isa VFKernel
        @test k.Q == [1.0, 2.0]
        @test size(k.R) == (2, 2)
        @test k.D == [0.5, 0.6]
    end

    @testset "chebykernel parsing and evaluate" begin
        d = Dict(
            "spectrum_bound" => 2.0,
            "approximations" => Any[
                Dict("coeffs" => [1.0, 2.0, 3.0]),
                Dict("coeffs" => [4.0, 5.0, 6.0]),
            ],
        )
        k = from_dict(ChebyKernel, d)
        @test k isa ChebyKernel
        @test size(k.C) == (3, 2)
        y = evaluate(k, [0.0, 1.0, 2.0])
        @test size(y) == (3, 2)
    end

    @testset "impulse defaults" begin
        L = path_laplacian(10)
        x = impulse(L, 3)
        @test ndims(x) == 1
        @test x[3] == 1.0
        @test sum(x) == 1.0

        X = impulse(L, 3, 2)
        @test size(X) == (10, 2)
        @test X[3, 1] == 1.0
        @test X[3, 2] == 1.0
    end

    @testset "spectral bound is positive" begin
        L = path_laplacian(20)
        ub = estimate_spectral_bound(L)
        @test ub > 0.01
    end

    @testset "resource registry basics" begin
        reg = SpectralGraphWavelet._ensure_registry()
        @test reg isa Dict
        @test haskey(reg, "MEXICAN_HAT")
        @test haskey(reg, "DELAY_TEXAS")
    end

    @testset "built-in graph and signals" begin
        L = resource("DELAY_TEXAS")
        @test issparse(L)
        @test size(L, 1) == size(L, 2)
        @test nnz(L) >= size(L, 1)

        S = resource("COORD_TEXAS")
        @test S isa AbstractArray
        @test ndims(S) == 2
        @test size(S, 2) in (2, 3)
    end

    @testset "mesh loading utilities" begin
        L = resource("MESH_BUNNY")
        xyz = resource("BUNNY_XYZ")
        @test issparse(L)
        @test size(L, 1) == size(xyz, 1)
        @test size(xyz, 2) == 3
    end

    @testset "list_graphs output" begin
        io = IOBuffer()
        list_graphs(io)
        out = String(take!(io))
        @test occursin("Graph Name", out)
        @test occursin("Vertices", out)
    end

    @testset "compatibility dll stubs" begin
        @test_throws ArgumentError get_cholmod_dll()
        @test_throws ArgumentError get_klu_dll()
    end
end
