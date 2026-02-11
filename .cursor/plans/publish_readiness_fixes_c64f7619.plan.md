---
name: Publish readiness fixes
overview: Fix the critical Cargo.toml errors and polish remaining items so the crate can be published to GitHub and crates.io.
todos:
  - id: fix-authors
    content: Fix Cargo.toml `authors` from string to array
    status: completed
  - id: fix-keywords
    content: Fix Cargo.toml `keywords` (remove spaces in 'long RNA reads')
    status: completed
  - id: fix-rustdoc
    content: Backtick-escape `<prefix>` and `<gene>` in 3 doc comments
    status: completed
  - id: fix-gitignore
    content: Replace Python .gitignore with Rust-appropriate version
    status: completed
  - id: add-lib-docs
    content: Add crate-level doc comment to src/lib.rs
    status: completed
  - id: fix-publish-md
    content: Update docs/PUBLISH.md placeholder to lrslab
    status: completed
  - id: update-changelog
    content: Add `flow` subcommand to CHANGELOG.md v0.1.0 entry
    status: completed
  - id: add-license-readme
    content: Re-add License section to README.md
    status: completed
  - id: exclude-toolchain
    content: Exclude rust-toolchain.toml from crate package
    status: completed
  - id: verify
    content: Run fmt, clippy, test, doc, and package checks
    status: completed
isProject: false
---

# Publish Readiness Fixes

## Critical (blocks all cargo commands)

### 1. Fix `authors` field in `Cargo.toml` (line 9)

Currently a bare string, which is invalid TOML for Cargo -- **every cargo command fails**:

```
authors = "runsheng"
```

Must be an array:

```toml
authors = ["runsheng"]
```

### 2. Fix `keywords` field in `Cargo.toml` (line 13)

`"long RNA reads"` contains spaces and exceeds single-keyword rules. crates.io rejects keywords with spaces. Replace with hyphenated or single-word alternatives:

```toml
keywords = ["bioinformatics", "genomics", "isoform", "nanopore", "long-reads"]
```

## Required (blocks clean publish)

### 3. Fix rustdoc warnings (3 files)

Angle-bracketed placeholders are parsed as HTML tags. Wrap each in backticks:

- [src/cli/preparedir.rs](src/cli/preparedir.rs) line 19: `<prefix>` (twice)
- [src/cli/flow.rs](src/cli/flow.rs) line 17: `<prefix>`
- [src/bin/clusterj_batch.rs](src/bin/clusterj_batch.rs) line 32: `<gene>`

### 4. Replace `.gitignore` with a Rust-appropriate version

The current [.gitignore](.gitignore) is a Python template (pycache, eggs, pip, etc.). It only has one Rust line (`/target/`). Replace contents with a proper Rust `.gitignore` that covers `/target/`, `Cargo.lock` exclusion for libraries (optional -- this is a binary crate so keeping it is correct), IDE files, `.DS_Store`, and the temp output patterns already present in the old version.

## Recommended (polish for a good first impression)

### 5. Add crate-level doc comment to `lib.rs`

[src/lib.rs](src/lib.rs) is currently 8 bare `pub mod` lines. Add a `//!` doc comment so docs.rs shows a meaningful landing page.

### 6. Update `docs/PUBLISH.md` placeholder

[docs/PUBLISH.md](docs/PUBLISH.md) line 32 still says `<your-user>`. Update to `lrslab` to match the Cargo.toml repository URL.

### 7. Update `CHANGELOG.md`

[CHANGELOG.md](CHANGELOG.md) doesn't mention the `flow` subcommand or `rust-toolchain.toml` pinning. Add these to the v0.1.0 entry.

### 8. Add License section back to `README.md`

The updated [README.md](README.md) lost its License section (previously present). Re-add it at the bottom.

### 9. Exclude `rust-toolchain.toml` from the crate package

The pinned `rust-toolchain.toml` (channel `1.90.0`) would force downstream `cargo install` users onto that exact toolchain. Add it to `exclude` in Cargo.toml, or add a `package.exclude` entry.

## Verification

After all fixes, run:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all --all-features
cargo doc --no-deps --all-features  # expect 0 warnings
cargo package --allow-dirty --no-verify
```

