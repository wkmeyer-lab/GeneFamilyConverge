# Installing the R packages (for the categorical phenotype tree step)

The **categorical phenotype tree** step
(`Src/Loc/scripts/make_categorical_phenotype_tree.R`) builds the CAFE-5 `-y`
(multi-λ) tree in R. It needs a small R stack:

| Package | Used for |
|---------|----------|
| `ape`    | Newick/`phylo` I/O and `reorder()` in the reconstruction |
| `castor` | `fit_mk()` — the Mk-model fit at the heart of the ancestral-state reconstruction |
| `expm`   | `expm()` — matrix exponential used in the pruning passes |
| `yaml`   | reading the pipeline `--config` file |

`Matrix` ships with R. **`phangorn` and `RERconverge` are NOT required**: the
reconstruction engine `getAncLiks()` is vendored in `Src/Reu/r/anc_recon.R` and
uses only `castor` + `expm` (with base R standing in for `phangorn::to.matrix`
and `Matrix::t`). This is byte-for-byte the lab's routine, just with a lighter
dependency footprint.

Like `tools/r8s`, the pipeline keeps **R and its packages out of the conda env**
(so the OrthoFinder/`convgeno` env stays lean). R is supplied by the system or,
on a cluster, by a module.

## Automated install (recommended)

No root needed — packages go to your personal R library.

```bash
module load R          # version is cluster-specific; skip on a workstation
Rscript tools/r-deps/install_r_packages.R
```

The script installs only what is missing (safe to re-run), auto-selects a
writable library (falling back to `~/R/library` when the system library is
read-only, as on most clusters), and verifies each package at the end.

Overrides:

```bash
# Install into a specific library:
Rscript tools/r-deps/install_r_packages.R "$HOME/R/library"

# Use a different CRAN mirror:
R_CRAN_MIRROR=https://cran.r-project.org Rscript tools/r-deps/install_r_packages.R
```

## conda-forge alternative

If you would rather manage R with conda, install it into a **separate** env (do
not add it to the `convgeno` env):

```bash
conda create -n convgeno-r -c conda-forge r-base r-ape r-castor r-expm r-yaml
conda activate convgeno-r
```

## Verify

```bash
Rscript -e 'for (p in c("ape","castor","expm","yaml")) cat(sprintf("%-8s %s\n", p, requireNamespace(p, quietly=TRUE)))'
```

Then the step runs as:

```bash
Rscript Src/Loc/scripts/make_categorical_phenotype_tree.R --config Src/Loc/configs/example_config.yaml
```

## Optional: exact RERconverge parity instead of the vendored engine

The vendored `getAncLiks()` already reproduces the lab's reconstruction. If you
specifically want to run the lab's installed `RERconverge` (heavier) or source
its `RERConvergeFunctions.R`, either install `RERconverge` or set
`phenotype_tree.rer_functions` in the config to the path of that file — it
overrides the vendored engine. Neither is needed for normal use.
