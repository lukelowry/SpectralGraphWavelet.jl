```@meta
CurrentModule = SpectralGraphWavelet
```

# SpectralGraphWavelet.jl

Sparse spectral graph wavelets and graph convolution in Julia.

## Installation

```julia
using Pkg
Pkg.add("SpectralGraphWavelet")
```

## Quick Start

```julia
using SpectralGraphWavelet

L = resource("DELAY_TEXAS")
x = impulse(L, 100)
conv = Convolve(L)
y = lowpass(conv, x, [0.1, 1.0, 10.0])
```

## Main Components

- `Convolve`, `DyConvolve`, `ChebyConvolve` for graph filtering
- `VFModel`, `ChebyModel` for fitting kernels
- `SGMA` for modal analysis
- `resource` and `resource_names` for bundled datasets/kernels

## API

```@autodocs
Modules = [SpectralGraphWavelet]
```
