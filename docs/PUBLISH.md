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
cargo package --allow-dirty --no-verify
cargo publish --dry-run
```

If `cargo publish --dry-run` fails due network/DNS in your environment, run it again from a network-enabled shell.

## 4) Publish to GitHub
```bash
git init
git add .
git commit -m "Initial release: trackcluster-rs v0.1.0"
git branch -M main
git remote add origin git@github.com:lrslab/trackcluster-rs.git
git push -u origin main
```

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
git tag v0.1.0
git push origin v0.1.0
```

This creates a GitHub Release with tarballs attached. Release notes are
generated automatically from commit history.
