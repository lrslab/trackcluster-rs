# Publishing Guide

This is the maintainer checklist for a GitHub binary release and the separate,
manual crates.io publication. Run it from a clean checkout of the commit that
will be tagged.

## 1. Freeze the release contents

- Choose the release version and update `Cargo.toml` and `Cargo.lock` together.
- Move the relevant `CHANGELOG.md` entries out of `Unreleased` into the new
  version section and review the migration/breaking-change notes.
- Confirm repository metadata, license files, README installation commands,
  supported targets, and minimum Rust version.
- Review the core release documentation:
  - `README.md`
  - `CHANGELOG.md`
  - `docs/CLI.md`
  - `docs/PIPELINE.md`
  - `docs/FORMATS.md`
  - `docs/MODIFICATION_VALIDATION.md`
  - `docs/INTERCHANGE.md`
  - `docs/behavior/`
- Run every README quickstart backed by the checked-in `examples/` inputs, and
  confirm the generated-BAM integration test covers the unbundled `bam2bigg`
  example. Public commands must not depend on private paths or `tests/fixtures`.
- Ensure all required docs and examples are tracked and that `git status` is
  clean. Do not commit build outputs, `.crate` files, SBOMs, checksums, or
  release tarballs.

Record the version that the release workflow will require:

```bash
version=$(cargo metadata --locked --no-deps --format-version 1 | jq -r '.packages[0].version')
test -n "$version"
```

## 2. Reproduce CI locally

The required local checks mirror `.github/workflows/ci.yml`:

```bash
cargo fmt --all -- --check
cargo clippy --locked --all-targets --all-features -- -D warnings
cargo test --locked --all --all-features
cargo bench --locked --all-features --no-run
RUSTDOCFLAGS="-D warnings" cargo doc --locked --no-deps --all-features
cargo package --locked
cargo audit
```

Also inspect the package contents and test crates.io packaging without
publishing:

```bash
cargo package --locked --list
cargo publish --locked --dry-run
```

The hosted CI additionally runs on Linux and macOS, rejects runtime Rust calls
to external processes, installs the generated `.crate`, smoke-tests both
binaries, and runs a one-thread flow. Both package and archive smoke tests
require a 40-hex Git commit. A package-path build must report a deterministic
`sha256:<64-hex>` source fingerprint, while the binary archive built directly
from the clean tagged checkout must report `source_fingerprint=clean`. This
also checks that Cargo's packaged VCS metadata preserves the tagged commit and
that editing an unpacked package cannot retain its cache identity. Do not tag
until the commit's complete CI run is green. If audit or dry-run needs network access, rerun it from a
network-enabled environment; do not waive it silently.

## 3. Create the GitHub release tag

The release workflow accepts tags matching `v*` and verifies that the tag is
exactly `v<Cargo package version>`. Push the release commit to `main`, wait for
CI, then create the tag from that exact commit:

```bash
test "$(git branch --show-current)" = "main"
test -z "$(git status --porcelain)"
test "v${version}" != "v"
git tag -a "v${version}" -m "Release v${version}"
git push origin "v${version}"
```

`.github/workflows/release.yml` builds and packages these targets:

- `x86_64-unknown-linux-musl`
- `aarch64-unknown-linux-gnu` (glibc 2.31+)
- `aarch64-apple-darwin`

Before any target build starts, the release workflow calls the reusable CI
workflow and requires its Linux/macOS formatting, Clippy, test, benchmark
compile, rustdoc, package-install smoke, and audit jobs to succeed.
The unpacked `x86_64-unknown-linux-musl` archive receives the end-to-end binary
smoke test; every target is checked for its required archive contents on its
build runner.

Each binary archive contains `trackcluster`, `clusterj_batch`,
`scripts/run_full_flow_rust.sh`, `README.md`, `CHANGELOG.md`, both license files,
the core CLI/pipeline/format/interchange and behavior documentation, and the
complete synthetic `examples/` directory. The workflow verifies this required
file set, rejects parent-directory entries, rejects README references to
repository-only fixtures/internal demos, and smoke-tests the unpacked Linux
archive, including the packaged full-flow script. The curated changelog section
is prepended to GitHub's generated commit/PR notes.

## 4. Verify release supply-chain artifacts

The GitHub release job also:

- generates `trackcluster-vX.Y.Z.cdx.json` in CycloneDX JSON format;
- writes `SHA256SUMS` for every binary tarball and the SBOM;
- creates GitHub build-provenance attestations covering the tarballs, SBOM, and
  checksum file;
- uploads all of those files to the GitHub Release.

Before announcing the release, verify that all three target tarballs, the
CycloneDX SBOM, `SHA256SUMS`, and provenance attestations are present. Download
the assets, validate the checksums, unpack at least one archive, and rerun
`trackcluster --version`, `trackcluster --help`, and
`clusterj_batch --version`.

## 5. Publish to crates.io manually

Pushing a Git tag does **not** publish the crate. crates.io publication is a
separate maintainer action and should use the same clean, tagged source after
the GitHub release succeeds:

```bash
git checkout "v${version}"
test -z "$(git status --porcelain)"
cargo publish --locked --dry-run
cargo login
cargo publish --locked
```

Confirm the new version and generated documentation on crates.io/docs.rs.
Publication cannot be replaced in place; if a release is wrong, yank it and
publish a new patch version rather than reusing the tag or version number.
