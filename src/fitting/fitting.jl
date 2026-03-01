module fitting

using LinearAlgebra
using Printf
using SparseArrays
using Statistics

import ..SpectralGraphWavelet: VFKernel, ChebyKernel, estimate_spectral_bound, _chebval_matrix

include("vfit.jl")
include("chebyshev.jl")

export VFModel, ChebyModel

end
