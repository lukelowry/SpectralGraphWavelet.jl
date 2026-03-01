using SpectralGraphWavelet
using Test

@testset "SpectralGraphWavelet.jl" begin
    include("test_helpers.jl")
    include("test_util.jl")
    include("test_functions.jl")
    include("test_vfit.jl")
    include("test_chebymodel.jl")
    include("test_chebyshev.jl")
    include("test_cholconv.jl")
    include("test_sgma.jl")
end
