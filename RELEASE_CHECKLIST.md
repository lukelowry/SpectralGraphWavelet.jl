# Release Checklist

## Versioning

- [ ] Bump `version` in `Project.toml` (next valid version after latest registry release).
- [ ] Keep semver meaning explicit in release notes.
- [ ] If breaking, include "breaking" in Registrator release notes.

## Compatibility

- [ ] Keep `[compat]` entries upper-bounded for all dependencies and `julia`.
- [ ] Ensure test-only deps in `[extras]` also have compat bounds.

## Artifacts and Data

- [ ] Verify every URL in `Artifacts.toml` is reachable from CI/clean environments.
- [ ] Confirm each artifact hash matches the hosted payload.
- [ ] Avoid relying on untracked local `data/` directories for release behavior.

## Validation

- [ ] `Pkg.test()` passes on minimum supported Julia (`1.10`) and latest stable.
- [ ] `Aqua.test_all` passes.
- [ ] Docs build passes (`julia --project=docs docs/make.jl`).

## Registry and Automation

- [ ] Ensure `TagBot` and `CompatHelper` secrets are configured (`DOCUMENTER_KEY`).
- [ ] Trigger Registrator from the intended commit/tag.
- [ ] Verify General PR checks pass before merge.
