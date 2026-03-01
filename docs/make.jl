using Documenter
using SpectralGraphWavelet

DocMeta.setdocmeta!(
    SpectralGraphWavelet,
    :DocTestSetup,
    :(using SpectralGraphWavelet),
    recursive=true,
)

makedocs(
    modules=[SpectralGraphWavelet],
    authors="luke-lowery",
    sitename="SpectralGraphWavelet.jl",
    format=Documenter.HTML(prettyurls=get(ENV, "CI", "false") == "true"),
    pages=[
        "Overview" => "index.md",
        "Usage" => "usage.md",
        "Library" => "library.md",
        "API Reference" => "api.md",
        "Related Projects" => "related.md",
    ],
)

deploydocs(
    repo="github.com/lukelowry/SpectralGraphWavelet.jl.git",
    devbranch="master",
)
