# Usage Guide

This page follows the same core flow as the Python `sgwt` usage guide, kept concise for Julia.

## Quickstart

```julia
using SpectralGraphWavelet

L = resource("DELAY_TEXAS")
signal = impulse(L, 600)
scales = [0.1, 1.0, 10.0]

conv = Convolve(L)
filtered = bandpass(conv, signal, scales)
```

### What Is Happening in This Example?

1. Load a built-in graph Laplacian (`DELAY_TEXAS`) and create an impulse signal.
2. Define spectral scales for multi-scale filtering.
3. Build a convolution engine (`Convolve`) and apply a spectral band-pass filter.
4. Receive one filtered output per requested scale.

## Underlying Graph

- Use bundled Laplacians from `resource("...")`.
- You can also load custom MAT/PLY data with `load_laplacian`, `load_ply_laplacian`, and `load_ply_xyz`.

## Input Signals

- `impulse(L, n)` creates a one-hot graph signal at vertex `n`.
- Signals should match graph size (`size(L, 1)`), with support for 1D and matrix forms.

## Kernel Functions

- Analytical kernels are available through `functions` and `convolve` helpers.
- Standard filters: `lowpass`, `bandpass`, `highpass`.

## Chebyshev Kernels

- Build polynomial approximations with `ChebyModel`.
- Apply approximated kernels via `ChebyConvolve`.

## Dynamic Graphs

```julia
using SpectralGraphWavelet

L = resource("DELAY_TEXAS")
poles = [10.0, 1.0, 0.1]
dconv = DyConvolve(L, poles)

x = impulse(L, 600)
y0 = lowpass(dconv, x)
addbranch!(dconv, 100, 200, 1.0)
y1 = lowpass(dconv, x)
```

## Notes

- Indices are 1-based in Julia (`impulse`, `addbranch!`, SGMA `bus` inputs).
- Scalar `scale` returns a single array; vector scales return multi-scale outputs.
- Python ground truth: <https://sgwt.readthedocs.io/>.
