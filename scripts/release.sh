#!/usr/bin/env bash
# Release builds for M-Tato (PROJECT.md milestone M9).
#
# What this script does on THIS machine (Linux x86-64):
#
#   1. `cargo build --release` for the two binaries, `mtato` and `mtato-gui`.
#   2. Copies them out of the workspace target directory into ./dist/linux-x86_64/.
#      The target directory is /tmp/mtato-target here (see .cargo/config.toml):
#      the root partition is 99% full, so build artefacts live on tmpfs and only
#      the finished binaries are copied onto disk.
#   3. Writes dist/VERSION and dist/SHA256SUMS.txt.
#
# What it deliberately does NOT do: cross-compile for Windows. There is no mingw
# toolchain on this machine and installing one needs root, so Windows x86-64
# binaries are produced by .github/workflows/release.yml instead. Pass --windows
# to try a cross build anyway; the script checks for the toolchain first and tells
# you exactly what is missing rather than failing half way through.
#
# Usage:
#   scripts/release.sh                 build for the host, stage into ./dist
#   scripts/release.sh --no-build      stage binaries that are already built
#   scripts/release.sh --check         report what would happen, change nothing
#   scripts/release.sh --windows       attempt a Windows cross build (needs mingw)
#   scripts/release.sh --target-dir DIR --dist DIR

set -euo pipefail

SCRIPT_PATH="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/$(basename -- "${BASH_SOURCE[0]}")"
REPO_ROOT="$(cd -- "$(dirname -- "$SCRIPT_PATH")/.." && pwd)"

BINARIES=(mtato mtato-gui)
CRATES=(mtato-cli mtato-gui)
HOST_TARGET="x86_64-unknown-linux-gnu"
WINDOWS_TARGET="x86_64-pc-windows-gnu"

DO_BUILD=1
DO_WINDOWS=0
CHECK_ONLY=0
TARGET_DIR=""
DIST_DIR="$REPO_ROOT/dist"

usage() {
    sed -n '2,30p' "$SCRIPT_PATH" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --no-build)   DO_BUILD=0 ;;
        --windows)    DO_WINDOWS=1 ;;
        --check)      CHECK_ONLY=1 ;;
        --target-dir) shift; TARGET_DIR="${1:?--target-dir needs a path}" ;;
        --dist)       shift; DIST_DIR="${1:?--dist needs a path}" ;;
        -h|--help)    usage 0 ;;
        *)            echo "ERROR: unknown argument: $1" >&2; usage 2 ;;
    esac
    shift
done

info() { printf '%s\n' "$*"; }
warn() { printf 'WARNING: %s\n' "$*" >&2; }
die()  { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# --------------------------------------------------------------------------
# Resolve the target directory. $CARGO_TARGET_DIR wins, then the value in
# .cargo/config.toml (which is what makes this /tmp/mtato-target here), then
# the cargo default ./target.
# --------------------------------------------------------------------------
resolve_target_dir() {
    if [ -n "$TARGET_DIR" ]; then
        printf '%s\n' "$TARGET_DIR"
        return
    fi
    if [ -n "${CARGO_TARGET_DIR:-}" ]; then
        printf '%s\n' "$CARGO_TARGET_DIR"
        return
    fi
    local config="$REPO_ROOT/.cargo/config.toml"
    if [ -f "$config" ]; then
        local value
        value="$(sed -n 's/^[[:space:]]*target-dir[[:space:]]*=[[:space:]]*"\(.*\)".*/\1/p' "$config" | head -n 1)"
        if [ -n "$value" ]; then
            printf '%s\n' "$value"
            return
        fi
    fi
    printf '%s\n' "$REPO_ROOT/target"
}

TARGET_DIR="$(resolve_target_dir)"
RELEASE_DIR="$TARGET_DIR/release"

version_of() {
    # Prefer the workspace version; `git describe` is only decoration and this
    # checkout may not be a git repository at all.
    local v
    v="$(sed -n 's/^version[[:space:]]*=[[:space:]]*"\(.*\)".*/\1/p' "$REPO_ROOT/Cargo.toml" | head -n 1)"
    [ -n "$v" ] || v="0.0.0"
    if command -v git >/dev/null 2>&1 && git -C "$REPO_ROOT" rev-parse --short HEAD >/dev/null 2>&1; then
        printf '%s+%s\n' "$v" "$(git -C "$REPO_ROOT" rev-parse --short HEAD)"
    else
        printf '%s\n' "$v"
    fi
}

VERSION="$(version_of)"

info "M-Tato release"
info "  repository:    $REPO_ROOT"
info "  version:       $VERSION"
info "  target dir:    $TARGET_DIR"
info "  dist dir:      $DIST_DIR"
info "  host target:   $HOST_TARGET"
info ""

# --------------------------------------------------------------------------
# Prerequisite checks: every one of them says how to fix the problem.
# --------------------------------------------------------------------------
if [ "$DO_BUILD" -eq 1 ]; then
    if ! command -v cargo >/dev/null 2>&1; then
        if [ -x "$HOME/.cargo/bin/cargo" ]; then
            warn "cargo is not on PATH but exists at \$HOME/.cargo/bin"
            export PATH="$HOME/.cargo/bin:$PATH"
        else
            die "cargo not found. Install Rust from https://rustup.rs, or run with --no-build
       to stage binaries that were built elsewhere."
        fi
    fi
fi

if [ "$DO_WINDOWS" -eq 1 ]; then
    MISSING=()
    if ! rustup target list --installed 2>/dev/null | grep -qx "$WINDOWS_TARGET"; then
        MISSING+=("rustup target add $WINDOWS_TARGET")
    fi
    if ! command -v x86_64-w64-mingw32-gcc >/dev/null 2>&1; then
        MISSING+=("a mingw-w64 cross linker (x86_64-w64-mingw32-gcc), which needs root to install")
    fi
    if [ "${#MISSING[@]}" -gt 0 ]; then
        info "Windows cross-compilation is not possible on this machine. Missing:"
        for m in "${MISSING[@]}"; do info "  - $m"; done
        info ""
        info "Windows x86-64 release binaries are built by CI instead:"
        info "  .github/workflows/release.yml"
        info "Push a tag (v*) or run the 'release' workflow from the Actions tab; the"
        info "binaries are attached to the run as artefacts and to the GitHub release."
        if [ "$CHECK_ONLY" -eq 1 ]; then
            exit 0
        fi
        die "refusing to pretend a Windows build happened; re-run without --windows
       for the Linux binaries, or use the GitHub Actions workflow for Windows."
    fi
fi

WORKFLOW="$REPO_ROOT/.github/workflows/release.yml"
if [ ! -f "$WORKFLOW" ]; then
    warn "no CI release workflow at $WORKFLOW - Windows binaries cannot be produced"
fi

if [ "$CHECK_ONLY" -eq 1 ]; then
    info "--check: nothing was built or copied."
    info "  would build:  cargo build --release --locked $(printf -- '--package %s ' "${CRATES[@]}")"
    info "  would copy:   $(printf '%s ' "${BINARIES[@]/#/$RELEASE_DIR/}")"
    info "  would write:  $DIST_DIR/$HOST_TARGET/, $DIST_DIR/VERSION, $DIST_DIR/SHA256SUMS.txt"
    exit 0
fi

# --------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------
if [ "$DO_BUILD" -eq 1 ]; then
    info "Building release binaries (this can take several minutes the first time)..."
    PKG_ARGS=()
    for c in "${CRATES[@]}"; do PKG_ARGS+=(--package "$c"); done
    ( cd "$REPO_ROOT" && cargo build --release --locked "${PKG_ARGS[@]}" )
    info ""
else
    info "--no-build: staging whatever is already in $RELEASE_DIR"
    info ""
fi

# --------------------------------------------------------------------------
# Stage
# --------------------------------------------------------------------------
stage() {
    local triple="$1" bindir="$2" exe_suffix="$3"
    local out="$DIST_DIR/$triple"
    mkdir -p "$out"
    local found=0
    for bin in "${BINARIES[@]}"; do
        local src="$bindir/$bin$exe_suffix"
        if [ ! -f "$src" ]; then
            warn "no $bin binary at $src"
            warn "     build it first: cargo build --release (target dir: $TARGET_DIR)"
            continue
        fi
        cp -f "$src" "$out/"
        chmod 0755 "$out/$(basename "$src")"
        info "  staged $(basename "$src") -> ${out#$REPO_ROOT/}"
        found=1
    done
    [ "$found" -eq 1 ] || die "nothing was staged for $triple: no binaries found in $bindir"
}

info "Staging into $DIST_DIR ..."
stage "$HOST_TARGET" "$RELEASE_DIR" ""
if [ "$DO_WINDOWS" -eq 1 ]; then
    ( cd "$REPO_ROOT" && cargo build --release --locked --target "$WINDOWS_TARGET" \
        --package mtato-cli --package mtato-gui )
    stage "$WINDOWS_TARGET" "$TARGET_DIR/$WINDOWS_TARGET/release" ".exe"
fi

printf '%s\n' "$VERSION" > "$DIST_DIR/VERSION"
info "  wrote VERSION ($VERSION)"

# Checksums, so a downloaded binary can be verified.
SUMS="$DIST_DIR/SHA256SUMS.txt"
: > "$SUMS"
if command -v sha256sum >/dev/null 2>&1; then
    ( cd "$DIST_DIR" && find . -type f ! -name SHA256SUMS.txt -print0 \
        | sort -z | xargs -0 sha256sum > "$SUMS" )
elif command -v shasum >/dev/null 2>&1; then
    ( cd "$DIST_DIR" && find . -type f ! -name SHA256SUMS.txt | sort \
        | xargs shasum -a 256 > "$SUMS" )
else
    warn "neither sha256sum nor shasum is available; $SUMS was left empty"
fi
info "  wrote SHA256SUMS.txt"

info ""
info "Done. Contents of $DIST_DIR:"
( cd "$DIST_DIR" && find . -type f | sort | sed 's/^/  /' )
info ""
info "Windows binaries: build them with .github/workflows/release.yml"
info "  (push a v* tag, or run the 'release' workflow from the Actions tab)"
