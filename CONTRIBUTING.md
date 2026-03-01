# Contributing

## Branching

1. Keep `master` (or `main`) as releasable.
2. Do all release-candidate work on a dedicated branch first.
3. Use short-lived topic branches for follow-up fixes.

Suggested branch names:

- `dev/prep-v0.1.1`
- `fix/<topic>`
- `docs/<topic>`

## Local Validation

Run package tests:

```powershell
julia --project=. -e "using Pkg; Pkg.test()"
```

Run quality checks:

```powershell
@'
using Pkg
Pkg.activate(; temp=true)
Pkg.develop(path=pwd())
Pkg.add("Aqua")
using Aqua, SpectralGraphWavelet
Aqua.test_all(SpectralGraphWavelet; persistent_tasks=false)
'@ | julia -
```

Build docs locally:

```powershell
julia --project=docs -e "using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()"
julia --project=docs docs/make.jl
```

## CI Workflows

- `CI.yml`: tests on supported Julia versions and OS targets.
- `Documentation.yml`: builds and deploys docs with Documenter.
- `CompatHelper.yml`: opens dependency compat PRs automatically.
- `TagBot.yml`: creates release tags and GitHub releases from registry merges.

## Release Process

Follow [`RELEASE_CHECKLIST.md`](RELEASE_CHECKLIST.md) before triggering Registrator.
