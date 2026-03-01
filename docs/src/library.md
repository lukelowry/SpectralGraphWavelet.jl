# Library

`SpectralGraphWavelet.jl` ships a bundled resource library for graph Laplacians, coordinates/signals, mesh assets, and analytical kernels.

## Discovering Resources

```julia
using SpectralGraphWavelet

names = resource_names()
list_graphs()
```

## Common Resource Keys

- Delays/Laplacians: `DELAY_TEXAS`, `DELAY_USA`, `DELAY_WECC`
- Impedance/Laplacians: `IMPEDANCE_TEXAS`, `IMPEDANCE_USA`
- Length/Laplacians: `LENGTH_TEXAS`, `LENGTH_USA`
- Coordinates/signals: `COORD_TEXAS`, `COORD_USA`
- Meshes: `MESH_BUNNY`, `MESH_HORSE`, `MESH_LBRAIN`, `MESH_RBRAIN`
- Mesh coordinates: `BUNNY_XYZ`, `HORSE_XYZ`, `LBRAIN_XYZ`, `RBRAIN_XYZ`
- Kernels: `MEXICAN_HAT`, `MODIFIED_MORLET`, `SHANNON`, `GAUSSIAN_WAV`

## Loading by Name

```julia
using SpectralGraphWavelet

L = resource("DELAY_TEXAS")
coords = resource("COORD_TEXAS")
kernel = resource("MEXICAN_HAT")
```

## Ground-Truth Alignment

The underlying library categories mirror the Python `sgwt` package library structure.

- Python project: <https://github.com/lukelowry/sgwt>
- Python docs library section: <https://sgwt.readthedocs.io/>
