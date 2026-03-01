module SpectralGraphWavelet

using LinearAlgebra
using MAT
using Printf
using SparseArrays
using SuiteSparse

export Convolve, DyConvolve, ChebyConvolve
export SGMA, ModeTable, NetworkAnalysisResult
export spectrum, analyze, analyze_many, find_peaks, find_modes, to_dict, to_array
export VFKernel, ChebyKernel, VFKern
export from_dict, evaluate
export VFModel, ChebyModel
export functions, fitting
export impulse, estimate_spectral_bound
export resource, resource_names
export convolve, convolve_multi, lowpass, bandpass, highpass
export addbranch, addbranch!
export load_laplacian, load_signal, load_ply_laplacian, load_ply_xyz, list_graphs
export get_cholmod_dll, get_klu_dll

include("util.jl")
include("functions/functions.jl")
include("fitting/fitting.jl")
include("cholconv.jl")
include("chebyconv.jl")
include("sgma.jl")

const VFModel = fitting.VFModel
const ChebyModel = fitting.ChebyModel

end
