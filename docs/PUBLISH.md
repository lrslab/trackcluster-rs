# Publishing Guide

This repository is prepared for both GitHub hosting and crates.io release.

## 1) Update repository metadata
Before first release, check the GitHub URL in `Cargo.toml`:
- `repository`
- `homepage`

## 2) Local validation (required)
```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all --all-features
```

## 3) Package checks
```bash
cargo package
cargo publish --dry-run
```

If `cargo publish --dry-run` fails due network/DNS in your environment, run it again from a network-enabled shell.

Do not commit build outputs or release tarballs to `main`; the release workflow uploads them to GitHub Releases as assets.

## 4) Publish to GitHub
Push a clean commit to `main` before tagging. Release artifacts should stay attached to GitHub Releases, not in the repository history.

## 5) Publish to crates.io
```bash
cargo login
cargo publish
```

## 6) Create GitHub release with binaries
The `.github/workflows/release.yml` workflow builds pre-built binaries for
Linux (x86_64 musl static, ARM64 glibc) and macOS (Apple Silicon)
automatically when a version tag is pushed:

```bash
git tag vX.Y.Z
git push origin vX.Y.Z
```

This creates a GitHub Release with tarballs attached. Release notes are
generated automatically from commit history.
