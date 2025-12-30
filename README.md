# SpectralGraphWavelet

[![Build Status](https://github.com/lukelowry/SpectralGraphWavelet.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/lukelowry/SpectralGraphWavelet.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Julia](https://img.shields.io/badge/julia-v1.10+-blue.svg)](https://julialang.org)

# Sparse GSP & SGWT Tools

A highly customizable, sparse-friendly SGWT/GSP module. This package provides tools to design, approximate, and implement a custom SGWT kernel for use over sparse networks.

Intended for Graph Signal Processing (GSP) of time-vertex signals over static and dynamic sparse graphs.

## Installation

```julia
using Pkg
Pkg.add("SpectralGraphWavelet")
```

## Usage

### Loading Data
The package provides utilities to load Laplacian matrices from `.mat` files.

```julia
using SpectralGraphWavelet

# Load graph topology (Laplacian)
# Expects a .mat file with key "A" or "L"
L = load_laplacian("path/to/graph.mat")
```

### Filtering Signals
The core of the package is the `DyConvolve` struct, which handles the Cholesky factorizations required for efficient filtering.

```julia
# Define scales/poles (poles are 1/s)
scale = 10.0
poles = [1.0 / scale] 

# Initialize convolution object
# This pre-computes Cholesky factorizations for the given poles
conv = DyConvolve(L, poles)

# Create an impulse signal at node 1
b = impulse(L, 1) 

# Apply filters
# These return a vector of results corresponding to the input poles
results_bp = bandpass(conv, b)
results_lp = lowpass(conv, b)
results_hp = highpass(conv, b)
```


### Dynamic Graph Updates
You can update the graph topology dynamically using `addbranch!`, which efficiently updates the internal factorizations using low-rank updates (via `SuiteSparse`).

```julia
# Update edge weight between node i and j by w
addbranch!(conv, i, j, w)
```


### Online Processing
For real-time applications, you can filter signals as they arrive and update the graph topology on the fly without full re-factorization.

```julia
# Initialize the engine once
conv = DyConvolve(L, poles)

# Simulate a data stream
for t in 1:1000
    # 1. Fetch the instantaneous signal vector at time t
    b_t = get_next_signal() 

    # 2. Filter the current time-step
    # The solver uses pre-computed factors, solve solve cost near linear in sparse system
    filtered_signal = lowpass(conv, b_t)

    # 3. Handle Dynamic Topology Changes
    # Efficiently update factors if the grid structure changes (e.g., a line closure)
    if t == 500
        # Add edge between i and j
        addbranch!(conv, i, j, weight) 
    end
end
```

## Theoretical Background

### Kernel Fitting

The kernel fitting representation is more generally a vector fitted function, a simple pole expansion of the form:
```math
g_a(\mathbf{\Lambda})\approx 
        d_aI + e_a\mathbf{\Lambda}
        + \sum_{q\in Q}\dfrac{r_{q,a}}{\mathbf{\Lambda}+qI} 
```

An iterative pole reallocation procedure is used to converge to a reduced order model. The convolution of some function $\mathbf{f}*g_a$ is computed using the Cholesky decomposition and memory efficient re-factors.

An example of an appropriate format of the rational expansion:

```json
{
    "nfunc": N,
    "d": [d0, d1, ..., dN],
    "npoles": M,
    "poles": [
        { "q": q0, "r":[r0, r1, ..., rN] },
        ...
    ]
}
```

### Analytical Filters

#### Low-Pass Spectral Graph Filter

The low-pass filter is *refinable*, as it is a self-similar rational function. The refinability makes it useful for signal smoothing across a range of spatial scales.

```math
\phi(\mathbf{\Lambda}) = \dfrac{I}{\mathbf{\Lambda}+I} 
```

#### High-Pass Spectral Graph Filter

The proposed high-pass filter acts as a container for variations over the graph below a given spatial scale.

```math
\mu(\mathbf{\Lambda}) = \dfrac{\mathbf{\Lambda}}{\mathbf{\Lambda}+I}
```

#### Band-Pass Spectral Graph Filter

A convenient closed-form wavelet generating kernel was found to be a useful kernel as an alternative to the vector-fitting procedure if a particular filter does not need to be designed. 

```math
\Psi(\mathbf{\Lambda}) = \dfrac{4\mathbf{\Lambda}}{(\mathbf{\Lambda}+I)^2} 
```

This filter qualifies as a wavelet generating kernel for the SGWT, since $\Psi(0)=0$ and the admissibility condition is satisfied. The admissibility constant of this band-pass filter is $C_f=8/3$.

### Cholesky Implementation

For the context that this was designed for, a direct solve approach is preferred to an iterative solver like `ARMA`. Time-varying graph signals must be as efficient as possible with memory, to ensure scalability to signals of large sparse networks. Especially if the process is online.

The `cholmod_solve2` and `updown` functions are the primary engine, in addition to other various design choices that accelerate graph convolutions.
