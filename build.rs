use std::ffi::OsString;
use std::fs;
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::process::Command;

use sha2::{Digest, Sha256};

fn json_string_field(contents: &str, key: &str) -> Option<String> {
    let key = format!("\"{key}\"");
    let after_key = contents.get(contents.find(&key)? + key.len()..)?;
    let after_colon = after_key.get(after_key.find(':')? + 1..)?.trim_start();
    let quoted = after_colon.strip_prefix('"')?;
    let value = quoted.get(..quoted.find('"')?)?;
    Some(value.to_owned())
}

struct PackagedVcs {
    commit: String,
}

fn packaged_vcs_info(manifest_dir: &Path) -> Option<PackagedVcs> {
    let contents = fs::read_to_string(manifest_dir.join(".cargo_vcs_info.json")).ok()?;
    let commit = json_string_field(&contents, "sha1")?;
    Some(PackagedVcs {
        commit: valid_commit(&commit)?,
    })
}

fn valid_commit(commit: &str) -> Option<String> {
    let commit = commit.trim();
    if commit.len() == 40 && commit.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        Some(commit.to_ascii_lowercase())
    } else {
        None
    }
}

fn git_output(manifest_dir: &Path, arguments: &[&str]) -> Option<Vec<u8>> {
    let output = Command::new("git")
        .arg("-C")
        .arg(manifest_dir)
        .env("GIT_OPTIONAL_LOCKS", "0")
        .args(arguments)
        .output()
        .ok()?;
    output.status.success().then_some(output.stdout)
}

fn git_commit(manifest_dir: &Path) -> Option<String> {
    let output = git_output(manifest_dir, &["rev-parse", "--verify", "HEAD"])?;
    valid_commit(std::str::from_utf8(&output).ok()?)
}

fn git_is_dirty(manifest_dir: &Path) -> Option<bool> {
    let output = git_output(
        manifest_dir,
        &[
            "status",
            "--porcelain=v1",
            "-z",
            "--untracked-files=all",
            "--",
            ".",
        ],
    )?;
    Some(!output.is_empty())
}

fn split_nul_paths(output: Vec<u8>) -> io::Result<Vec<PathBuf>> {
    let mut paths = output
        .split(|byte| *byte == 0)
        .filter(|path| !path.is_empty())
        .map(bytes_to_path)
        .collect::<io::Result<Vec<_>>>()?;
    paths.sort_by_key(|path| path_bytes(path));
    paths.dedup();
    Ok(paths)
}

#[cfg(unix)]
fn bytes_to_path(bytes: &[u8]) -> io::Result<PathBuf> {
    use std::os::unix::ffi::OsStringExt;

    Ok(PathBuf::from(OsString::from_vec(bytes.to_vec())))
}

#[cfg(not(unix))]
fn bytes_to_path(bytes: &[u8]) -> io::Result<PathBuf> {
    let path = String::from_utf8(bytes.to_vec())
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    Ok(PathBuf::from(path))
}

#[cfg(unix)]
fn path_bytes(path: &Path) -> Vec<u8> {
    use std::os::unix::ffi::OsStrExt;

    path.as_os_str().as_bytes().to_vec()
}

#[cfg(not(unix))]
fn path_bytes(path: &Path) -> Vec<u8> {
    path.to_string_lossy().replace('\\', "/").into_bytes()
}

#[cfg(unix)]
fn is_executable(metadata: &fs::Metadata) -> bool {
    use std::os::unix::fs::PermissionsExt;

    metadata.permissions().mode() & 0o111 != 0
}

#[cfg(not(unix))]
fn is_executable(_metadata: &fs::Metadata) -> bool {
    false
}

fn git_source_paths(manifest_dir: &Path) -> Option<Vec<PathBuf>> {
    let output = git_output(
        manifest_dir,
        &[
            "ls-files",
            "--cached",
            "--others",
            "--exclude-standard",
            "-z",
            "--",
            ".",
        ],
    )?;
    split_nul_paths(output).ok()
}

fn packaged_source_paths(manifest_dir: &Path) -> io::Result<Vec<PathBuf>> {
    fn visit(root: &Path, relative: &Path, paths: &mut Vec<PathBuf>) -> io::Result<()> {
        let directory = root.join(relative);
        let mut entries = fs::read_dir(&directory)?.collect::<Result<Vec<_>, _>>()?;
        entries.sort_by_key(|entry| {
            let name = entry.file_name();
            path_bytes(Path::new(&name))
        });
        for entry in entries {
            let name = entry.file_name();
            if relative.as_os_str().is_empty()
                && matches!(
                    name.to_str(),
                    Some(".git" | "target" | ".cargo-ok" | ".cargo-checksum.json")
                )
            {
                continue;
            }
            let child = relative.join(name);
            let metadata = fs::symlink_metadata(entry.path())?;
            if metadata.is_dir() && !metadata.file_type().is_symlink() {
                visit(root, &child, paths)?;
            } else {
                paths.push(child);
            }
        }
        Ok(())
    }

    let mut paths = Vec::new();
    visit(manifest_dir, Path::new(""), &mut paths)?;
    paths.sort_by_key(|path| path_bytes(path));
    Ok(paths)
}

fn source_fingerprint(manifest_dir: &Path, paths: &[PathBuf]) -> io::Result<String> {
    let mut hasher = Sha256::new();
    hasher.update(b"trackcluster-source-fingerprint-v1\0");
    let mut buffer = [0u8; 64 * 1024];

    for relative in paths {
        let encoded_path = path_bytes(relative);
        hasher.update((encoded_path.len() as u64).to_le_bytes());
        hasher.update(&encoded_path);
        let path = manifest_dir.join(relative);
        match fs::symlink_metadata(&path) {
            Ok(metadata) if metadata.file_type().is_symlink() => {
                hasher.update(b"symlink\0");
                let target = fs::read_link(&path)?;
                let target = path_bytes(&target);
                hasher.update((target.len() as u64).to_le_bytes());
                hasher.update(target);
            }
            Ok(metadata) if metadata.is_file() => {
                hasher.update(b"file\0");
                hasher.update([u8::from(is_executable(&metadata))]);
                hasher.update(metadata.len().to_le_bytes());
                let mut file = fs::File::open(&path)?;
                loop {
                    let read = file.read(&mut buffer)?;
                    if read == 0 {
                        break;
                    }
                    hasher.update(&buffer[..read]);
                }
            }
            Ok(_) => hasher.update(b"other\0"),
            Err(error) if error.kind() == io::ErrorKind::NotFound => {
                // Deleted tracked paths remain in `git ls-files` and therefore
                // deliberately affect a dirty source snapshot.
                hasher.update(b"missing\0");
            }
            Err(error) => return Err(error),
        }
    }

    Ok(format!("sha256:{:x}", hasher.finalize()))
}

fn emit_rerun_paths(manifest_dir: &Path, paths: &[PathBuf], watch_git: bool) {
    for relative in paths {
        if !path_bytes(relative).contains(&b'\n')
            && fs::symlink_metadata(manifest_dir.join(relative)).is_ok()
        {
            println!(
                "cargo:rerun-if-changed={}",
                manifest_dir.join(relative).display()
            );
        }
    }
    if !watch_git {
        return;
    }
    for git_path in ["HEAD", "index", "packed-refs", "refs"] {
        let Some(output) = git_output(manifest_dir, &["rev-parse", "--git-path", git_path]) else {
            continue;
        };
        let Ok(path) = std::str::from_utf8(&output) else {
            continue;
        };
        let path = PathBuf::from(path.trim());
        let path = if path.is_absolute() {
            path
        } else {
            manifest_dir.join(path)
        };
        println!("cargo:rerun-if-changed={}", path.display());
    }
}

fn nonempty_override(name: &str) -> Option<String> {
    let value = std::env::var(name).ok()?;
    let value = value.trim();
    (!value.is_empty() && !value.contains(['\r', '\n'])).then(|| value.to_owned())
}

fn main() {
    println!("cargo:rerun-if-env-changed=TRACKCLUSTER_GIT_COMMIT");
    println!("cargo:rerun-if-env-changed=TRACKCLUSTER_SOURCE_FINGERPRINT");

    let manifest_dir = std::env::var_os("CARGO_MANIFEST_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));
    let packaged_vcs = packaged_vcs_info(&manifest_dir);
    let source_paths = if packaged_vcs.is_some() {
        packaged_source_paths(&manifest_dir)
    } else {
        git_source_paths(&manifest_dir)
            .map(Ok)
            .unwrap_or_else(|| packaged_source_paths(&manifest_dir))
    }
    .unwrap_or_else(|error| panic!("cannot enumerate source files for build identity: {error}"));
    emit_rerun_paths(&manifest_dir, &source_paths, packaged_vcs.is_none());

    let commit = nonempty_override("TRACKCLUSTER_GIT_COMMIT")
        .or_else(|| packaged_vcs.as_ref().map(|vcs| vcs.commit.clone()))
        .or_else(|| git_commit(&manifest_dir))
        .unwrap_or_else(|| "unknown".to_owned());
    let fingerprint = nonempty_override("TRACKCLUSTER_SOURCE_FINGERPRINT").unwrap_or_else(|| {
        if packaged_vcs.is_some() {
            // A local `.crate` archive has VCS metadata but no immutable
            // per-file checksum manifest. Hash the actual packaged tree on
            // every build so edits made after unpacking cannot inherit the
            // official package's cache identity. Registry-generated checksum
            // metadata is excluded from `packaged_source_paths` above.
            source_fingerprint(&manifest_dir, &source_paths)
                .unwrap_or_else(|error| panic!("cannot fingerprint build source: {error}"))
        } else if git_is_dirty(&manifest_dir).is_some_and(|dirty| !dirty) {
            "clean".to_owned()
        } else {
            source_fingerprint(&manifest_dir, &source_paths)
                .unwrap_or_else(|error| panic!("cannot fingerprint build source: {error}"))
        }
    });
    println!("cargo:rustc-env=TRACKCLUSTER_GIT_COMMIT={commit}");
    println!("cargo:rustc-env=TRACKCLUSTER_SOURCE_FINGERPRINT={fingerprint}");
}
