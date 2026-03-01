```@meta
CurrentModule = SpectralGraphWavelet
```

# API Reference

This reference follows the same high-level grouping as the Python `sgwt` API.

## Graph Convolution

- Engines: `Convolve`, `DyConvolve`, `ChebyConvolve`
- Core operations: `convolve`, `lowpass`, `bandpass`, `highpass`
- Dynamic updates: `addbranch`, `addbranch!`

## Applications

- Modal analysis: `SGMA`, `ModeTable`
- Analysis functions: `spectrum`, `analyze`, `analyze_many`, `find_peaks`, `find_modes`

## Utility and Resources

- Signal and spectrum helpers: `impulse`, `estimate_spectral_bound`
- Resource loading: `resource`, `resource_names`, `list_graphs`
- File loaders: `load_laplacian`, `load_signal`, `load_ply_laplacian`, `load_ply_xyz`

## Kernel Fitting and Functions

```@docs
SpectralGraphWavelet.functions.gaussian_wavelet
SpectralGraphWavelet.functions.lowpass
SpectralGraphWavelet.functions.bandpass
SpectralGraphWavelet.functions.highpass
```

```@autodocs
Modules = [SpectralGraphWavelet]
Pages = ["util.jl", "fitting/fitting.jl", "functions/functions.jl", "chebyconv.jl", "cholconv.jl", "sgma.jl"]
```
