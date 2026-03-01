# SpectralGraphWavelet.jl

Julia port of the `sgwt` package for sparse spectral graph wavelets and graph convolution, using native Julia `SuiteSparse` (no DLL wrappers required).

## Scope

This package mirrors the Python `sgwt` module surface for:

- Static and dynamic sparse graph convolution (`Convolve`, `DyConvolve`)
- Chebyshev convolution (`ChebyConvolve`)
- Vector-fitting and Chebyshev fitting models (`VFModel`, `ChebyModel`)
- SGMA modal analysis (`SGMA`, `ModeTable`)
- Analytical kernels (`SpectralGraphWavelet.functions`)
- Bundled library assets (artifact-backed, loaded via `resource("NAME")`)

## Install

```julia
using Pkg
Pkg.add("SpectralGraphWavelet")
```

## Documentation

This repository uses the standard Julia `Documenter.jl` workflow.

- Source: `docs/src/index.md`
- Builder: `docs/make.jl`
- CI workflow: `.github/workflows/Documentation.yml`

Build docs locally:

```julia
using Pkg
Pkg.activate("docs")
Pkg.instantiate()
include("docs/make.jl")
```

## Quick Start

```julia
using SpectralGraphWavelet

L = resource("DELAY_TEXAS")
x = impulse(L, 100)

conv = Convolve(L)
y_lp = lowpass(conv, x, [0.1, 1.0, 10.0])
y_bp = bandpass(conv, x, 1.0)   # scalar scale returns a single array
```

## Dynamic Topology

```julia
poles = [10.0, 1.0, 0.1]
dconv = DyConvolve(resource("DELAY_TEXAS"), poles)

y0 = lowpass(dconv, impulse(resource("DELAY_TEXAS"), 100))
addbranch(dconv, 100, 200, 1.0)
y1 = lowpass(dconv, impulse(resource("DELAY_TEXAS"), 100))
```

## Fitting + Convolution

```julia
lam = collect(range(0.01, 100.0, length=500))
g = @. 1 / (lam + 1)

K = VFModel.kernel(lam, g, 6)
conv = Convolve(resource("DELAY_TEXAS"))
Y = convolve(conv, impulse(resource("DELAY_TEXAS"), 100), K)
```

## SGMA

```julia
L = resource("DELAY_HAWAII")
scales = exp10.(range(-1, 1, length=10))
freqs = collect(range(0.1, 2.0, length=20))
engine = SGMA(L, scales, freqs)
```

## Resources

Resources are discovered from the bundled artifact and loaded on demand with `resource("NAME")`, for example:

- `DELAY_TEXAS`, `DELAY_USA`, `DELAY_WECC`
- `IMPEDANCE_TEXAS`, `LENGTH_USA`
- `COORD_TEXAS`, `COORD_USA`
- `MESH_BUNNY`, `BUNNY_XYZ`
- `MEXICAN_HAT`, `MODIFIED_MORLET`, `SHANNON`

Use `resource_names()` to list discovered resource keys and `list_graphs()` to print graph assets.

## Python Compatibility Notes

- No `dll` wrappers are shipped in Julia; `get_cholmod_dll()` and `get_klu_dll()` are compatibility stubs and intentionally throw.
- Indices are Julia-native 1-based in all APIs (`impulse`, `addbranch`, SGMA `bus`).
- `Convolve`/`DyConvolve` retain the Python-style signal handling semantics: complex split/recombine, scalar scales, and 1D signal squeeze behavior.

## Development

- Contribution guide: `CONTRIBUTING.md`
- Release checklist: `RELEASE_CHECKLIST.md`
