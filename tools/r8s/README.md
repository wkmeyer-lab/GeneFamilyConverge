# Installing r8s (for the ultrametric step)

The pipeline uses **r8s 1.81** to turn OrthoFinder's species tree into the
rooted, ultrametric time tree CAFE-5 requires (see
`Src/Loc/scripts/make_tree_ultrametric.py`).

r8s is **not on conda/bioconda** — it is distributed only as 2004-era source
(SourceForge / ginger.ucdavis.edu), and its bundled binary is macOS-only, so on
Linux you must compile it. The source also predates modern GCC defaults, so a
plain `make` fails. `install_r8s.sh` in this directory applies the needed fixes
and builds a working binary; this README is the manual fallback for clusters
where the automated build needs adjusting.

## Automated install (recommended)

```bash
# 1. Load a modern GCC toolchain (version is cluster-specific).
module spider gcc 2>/dev/null || module avail gcc
module load gcc/12.4.0

# 2. Build + install to $HOME/.local/bin/r8s.
chmod +x tools/r8s/install_r8s.sh
./tools/r8s/install_r8s.sh
```

Custom locations via environment variables:

```bash
R8S_INSTALL_ROOT="$HOME/apps" R8S_BIN_DIR="$HOME/bin" ./tools/r8s/install_r8s.sh
```

The installer downloads the official source, removes the bundled macOS binary
(keeping the required `TN/` directory), uses the lowercase `makefile` (the
archive's `Makefile.linux` wrongly requires the absent `continuousML.h`), adds
`-fcommon` (GCC ≥10) and `-std=legacy -fallow-argument-mismatch` (legacy
Fortran), embeds the `libgfortran` directory as an ELF rpath, and verifies the
result is a Linux ELF executable reporting `r8s version 1.81`.

## Point the pipeline at it

The pipeline resolves r8s from `Src/Loc/configs/tool_paths.yaml`
(`r8s.command`, default `r8s`). Choose one:

- **On PATH:** add `$HOME/.local/bin` to `PATH` (the installer prints the line;
  it also belongs in your SLURM job scripts), then the default `r8s.command:
  r8s` just works.
- **Absolute path:** set it in the gitignored `Src/Loc/configs/tool_paths.local.yaml`:

  ```yaml
  r8s:
    command: /home/USER/.local/bin/r8s
  ```

Or pass `--r8s-path /home/USER/.local/bin/r8s` to `make_tree_ultrametric.py`.

## Verify

```bash
command -v r8s
r8s -v -b                       # -> r8s version 1.81
file "$(command -v r8s)"        # -> ELF 64-bit ... executable  (NOT Mach-O)

# rpath check: must still print the version with LD_LIBRARY_PATH cleared,
# otherwise the binary can't find libgfortran on a fresh compute node.
env -u LD_LIBRARY_PATH r8s -v -b
```

If the rpath check fails with a missing `libgfortran.so`, either `module load`
the same GCC in every SLURM script, or rebuild with that module loaded.

## Troubleshooting

- **`gcc`/`gfortran` not found** — `module spider gcc` then load a version.
- **Existing binary reports `Mach-O`** — that's the archive's macOS binary;
  delete only the `r8s` file (never `TN/`) and rebuild.
- **Build warnings** — expected for this old source; success is exit code `0`
  plus `r8s version 1.81` and an `ELF` file.
- **Non-GNU toolchains** (Intel/NVIDIA/AOCC), non-glibc systems, or blocked
  SourceForge access — the automated flags target GCC/gfortran; these cases may
  need manual compiler-specific changes.

## Alternative: container

If your cluster supports Apptainer/Singularity, a container with r8s avoids the
per-cluster build entirely; point `r8s.command` at a wrapper that invokes
`apptainer exec <image> r8s`. The source build above is the default because it
needs no container runtime.
