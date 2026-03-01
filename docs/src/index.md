```@meta
CurrentModule = SpectralGraphWavelet
```

# SpectralGraphWavelet.jl

Sparse graph signal processing and spectral graph wavelets in Julia.

This package is the Julia implementation aligned to the Python [`sgwt`](https://github.com/lukelowry/sgwt) project.

For conceptual parity and examples, treat the Python docs as ground truth:

- `sgwt` docs: <https://sgwt.readthedocs.io/>
- `sgwt` repository: <https://github.com/lukelowry/sgwt>
- Julia repository: <https://github.com/lukelowry/SpectralGraphWavelet.jl>
- Research website: <https://lukelowry.github.io/>

## Installation

```julia
using Pkg
Pkg.add("SpectralGraphWavelet")
```

## Quickstart

```julia
using SpectralGraphWavelet

L = resource("DELAY_TEXAS")
x = impulse(L, 100)
conv = Convolve(L) # static graph convolution engine
y = bandpass(conv, x, [0.1, 1.0, 10.0])
```

## Package Scope

- Static graph filtering: `Convolve`
- Dynamic topology filtering: `DyConvolve`
- Chebyshev approximation filtering: `ChebyConvolve`
- Kernel fitting: `VFModel`, `ChebyModel`
- Modal analysis: `SGMA`
- Built-in resource library: `resource`, `resource_names`, `list_graphs`

See:

- [Usage](usage.md)
- [Library](library.md)
- [API Reference](api.md)
