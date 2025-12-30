module SpectralGraphWavelet

using MAT
using SparseArrays
using LinearAlgebra
using SuiteSparse
using Statistics

export load_laplacian, load_signal
export VFKern, DyConvolve, impulse, convolve, lowpass, bandpass, highpass, addbranch!

include("io.jl")
include("dynamic.jl")

end
