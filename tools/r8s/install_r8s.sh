#!/usr/bin/env bash
# Install r8s 1.81 from source into a user-owned directory on a modern GNU/Linux
# system. No root privileges are required.
#
# r8s is NOT available on conda/bioconda and its 2004-era source needs a few
# fixes to compile with modern GCC/gfortran (see tools/r8s/README.md). This
# installer applies those fixes, builds an ELF binary, embeds the libgfortran
# runtime path (rpath) so it survives across SLURM sessions, and installs to
# $HOME/.local/bin/r8s.
#
# Usage (load a modern GCC module first):
#   module load gcc/12.4.0        # version is cluster-specific
#   ./tools/r8s/install_r8s.sh
#
# Then either put $HOME/.local/bin on PATH, or set the absolute path in
# Src/Loc/configs/tool_paths.local.yaml (r8s.command).
#
# Tested on a Lehigh SLURM cluster with an environment-modules GCC toolchain.

set -Eeuo pipefail

R8S_VERSION="1.81"
INSTALL_ROOT="${R8S_INSTALL_ROOT:-$HOME/software}"
BIN_DIR="${R8S_BIN_DIR:-$HOME/.local/bin}"
ARCHIVE="${INSTALL_ROOT}/r8s${R8S_VERSION}.tar.gz"
SOURCE_DIR="${INSTALL_ROOT}/r8s${R8S_VERSION}"
DOWNLOAD_URL="https://sourceforge.net/projects/r8s/files/r8s${R8S_VERSION}.tar.gz/download"

log() {
    printf '[r8s installer] %s\n' "$*"
}

fail() {
    printf '[r8s installer] ERROR: %s\n' "$*" >&2
    exit 1
}

[[ "$(uname -s)" == "Linux" ]] || fail "This installer is for Linux."

# Respect compiler variables set by a module or compiler environment.
CC_BIN="${CC:-gcc}"
FC_BIN="${FC:-gfortran}"

command -v "$CC_BIN" >/dev/null 2>&1 || \
    fail "C compiler '$CC_BIN' was not found. Load a GCC module or set CC."
command -v "$FC_BIN" >/dev/null 2>&1 || \
    fail "Fortran compiler '$FC_BIN' was not found. Load a GCC/gfortran module or set FC."
command -v make >/dev/null 2>&1 || fail "GNU make was not found."
command -v tar >/dev/null 2>&1 || fail "tar was not found."
command -v gzip >/dev/null 2>&1 || fail "gzip was not found."
command -v install >/dev/null 2>&1 || fail "the install command was not found."

if command -v curl >/dev/null 2>&1; then
    DOWNLOADER=(curl --fail --location --retry 3 --output "$ARCHIVE" "$DOWNLOAD_URL")
elif command -v wget >/dev/null 2>&1; then
    DOWNLOADER=(wget --tries=3 --output-document="$ARCHIVE" "$DOWNLOAD_URL")
else
    fail "Neither curl nor wget was found."
fi

log "Using C compiler: $(command -v "$CC_BIN")"
log "Using Fortran compiler: $(command -v "$FC_BIN")"
log "Installing under: $INSTALL_ROOT"

mkdir -p "$INSTALL_ROOT" "$BIN_DIR"

log "Downloading r8s ${R8S_VERSION} from SourceForge"
"${DOWNLOADER[@]}"
[[ -s "$ARCHIVE" ]] || fail "The downloaded archive is empty."

gzip -t "$ARCHIVE" || fail "The downloaded archive is not a valid gzip file."

log "Extracting source archive"
rm -rf "$SOURCE_DIR"
tar -xzf "$ARCHIVE" -C "$INSTALL_ROOT"
SRC="${SOURCE_DIR}/src"
[[ -d "$SRC" ]] || fail "Expected source directory was not found: $SRC"
cd "$SRC"

# The archive includes a macOS executable named r8s. Remove only that file.
# TN is a required source directory and must not be removed.
rm -f ./r8s ./*.o
[[ -d TN ]] || fail "Required TN source directory is missing."
[[ -f makefile ]] || fail "The lowercase source makefile is missing."

cp makefile makefile.original

# The lowercase makefile excludes the unavailable continuousML source. The
# archive's Makefile.linux incorrectly requires continuousML.h, so do not use it.

# Remove obsolete hard-coded /usr/include dependencies. Modern compilers locate
# standard headers through their configured include search paths.
sed -E -i 's#/usr/include/[^[:space:]]+##g' makefile

# r8s predates GCC 10's default change to -fno-common.
if grep -qE '^CFLAGS[[:space:]]*=' makefile; then
    sed -E -i 's|^CFLAGS[[:space:]]*=.*$|CFLAGS = -O2 -fcommon|' makefile
else
    printf '\nCFLAGS = -O2 -fcommon\n' >> makefile
fi

# Permit legacy Fortran constructs and historical argument mismatches.
if grep -qE '^FFLAGS[[:space:]]*=' makefile; then
    sed -E -i \
        's|^FFLAGS[[:space:]]*=.*$|FFLAGS = -O2 -std=legacy -fallow-argument-mismatch|' \
        makefile
else
    printf '\nFFLAGS = -O2 -std=legacy -fallow-argument-mismatch\n' >> makefile
fi

# Locate libgfortran from the selected compiler. Embed a runtime search path so
# the installed binary can locate the same compiler runtime in later sessions.
GFORTRAN_LIBRARY="$($FC_BIN --print-file-name=libgfortran.so)"
[[ "$GFORTRAN_LIBRARY" = /* && -e "$GFORTRAN_LIBRARY" ]] || \
    fail "Could not locate libgfortran.so using '$FC_BIN --print-file-name'."
GFORTRAN_LIB_DIR="$(dirname "$GFORTRAN_LIBRARY")"

if grep -qE '^LPATH[[:space:]]*=' makefile; then
    sed -E -i \
        "s|^LPATH[[:space:]]*=.*$|LPATH = -L${GFORTRAN_LIB_DIR} -Wl,-rpath,${GFORTRAN_LIB_DIR}|" \
        makefile
else
    printf '\nLPATH = -L%s -Wl,-rpath,%s\n' \
        "$GFORTRAN_LIB_DIR" "$GFORTRAN_LIB_DIR" >> makefile
fi

log "Building r8s"
make -f makefile \
    CC="$CC_BIN" \
    FC="$FC_BIN" \
    F77="$FC_BIN"

[[ -x ./r8s ]] || fail "Build finished without creating an executable r8s file."

VERSION_OUTPUT="$(./r8s -v -b 2>&1 || true)"
printf '%s\n' "$VERSION_OUTPUT"
grep -q "r8s version ${R8S_VERSION}" <<<"$VERSION_OUTPUT" || \
    fail "The compiled program did not report r8s version ${R8S_VERSION}."

if command -v file >/dev/null 2>&1; then
    FILE_OUTPUT="$(file ./r8s)"
    printf '%s\n' "$FILE_OUTPUT"
    grep -q 'ELF' <<<"$FILE_OUTPUT" || fail "The compiled file is not a Linux ELF executable."
fi

if command -v ldd >/dev/null 2>&1; then
    if ldd ./r8s 2>&1 | grep -q 'not found'; then
        ldd ./r8s >&2 || true
        fail "The compiled executable has unresolved shared-library dependencies."
    fi
fi

log "Installing executable to $BIN_DIR/r8s"
install -m 0755 ./r8s "$BIN_DIR/r8s"

log "Installation succeeded"
printf '\nExecutable: %s\n' "$BIN_DIR/r8s"
printf 'Version:    %s\n' "$("$BIN_DIR/r8s" -v -b 2>&1 | head -n 1)"
printf '\nAdd this line to your shell startup file if %s is not already on PATH:\n' "$BIN_DIR"
printf 'export PATH="%s:$PATH"\n' "$BIN_DIR"
printf '\nThen point the pipeline at it (if not on PATH), in\n'
printf 'Src/Loc/configs/tool_paths.local.yaml:\n'
printf 'r8s:\n  command: %s/r8s\n' "$BIN_DIR"
