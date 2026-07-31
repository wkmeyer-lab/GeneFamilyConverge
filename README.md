# GeneFamilyConverge

**A pipeline for characterizing gene-family copy number across species and
associating copy-number variation with convergent phenotypes.**

GeneFamilyConverge takes a set of proteomes and walks them through the standard
comparative-genomics chain — clean the proteomes, infer gene families, build a
time-calibrated species tree, model gene-family expansion/contraction, and
associate copy-number change with convergent traits. Its distinguishing feature
is a **multi-node, SLURM-orchestrated OrthoFinder workflow** that spreads the
expensive all-vs-all protein search across a whole cluster while producing
results identical to a single-node OrthoFinder run.

The pipeline is command-line driven through the `convgeno` tool plus a small set
of stage scripts. This README tells you how to install it, what each stage does,
and how to run them.

---

## Contents

- [The pipeline at a glance](#the-pipeline-at-a-glance)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Usage](#usage)
  - [1. Filter isoforms](#1-filter-isoforms)
  - [2. Set up your run (`convgeno init`)](#2-set-up-your-run-convgeno-init)
  - [3. Run OrthoFinder and build the time-calibrated tree](#3-run-orthofinder-and-build-the-time-calibrated-tree)
  - [4. Model gene-family turnover with CAFE-5](#4-model-gene-family-turnover-with-cafe-5)
  - [5. Associate copy number with convergent traits](#5-associate-copy-number-with-convergent-traits)
  - [Calibrating DIAMOND throughput](#calibrating-diamond-throughput)
- [Inputs and outputs](#inputs-and-outputs)
- [Configuration](#configuration)
- [Command reference](#command-reference)
- [Documentation](#documentation)
- [How to cite](#how-to-cite)
- [Acknowledgements](#acknowledgements)

---

## The pipeline at a glance

```
   proteomes                one FASTA per species
       │
       ▼
 (1) filter isoforms        longest protein per gene            convgeno / filter_isoforms.py
       │
       ▼
 (2) OrthoFinder            gene families (orthogroups)         convgeno orthofinder
       │  + species tree
       ▼
 (2b) r8s                   ultrametric, time-calibrated        (runs automatically with OrthoFinder)
       │                    species tree
       ▼
 (3) CAFE-5                 gene-family expansion / contraction
       │
       ▼
 (4) R analysis             copy-number change ↔ convergent      RERconverge-based association
                            phenotype
```

The scientific question is whether **shifts in gene-family copy number** line up
with **independently evolved (convergent) phenotypes** — for example, whether a
gene family repeatedly expands in lineages that share a diet or habitat.
OrthoFinder defines the gene families, r8s scales the species tree to time,
CAFE-5 models where families expand and contract, and a downstream R step links
those changes to convergent traits.

---

## Requirements

**Core pipeline (proteome filtering and OrthoFinder):**

- A **SLURM** HPC cluster with `sinfo`, `scontrol`, and `sbatch` on the login
  node (`sacctmgr` is used opportunistically if present). The OrthoFinder
  workflow submits cluster jobs; the isoform filter can run anywhere.
- **conda** (miniforge or miniconda).
- **Python 3.10–3.12**.
- Bioinformatics tools, all installed for you by the conda environment:
  **OrthoFinder 2.5.5**, **DIAMOND 2.1.9**, **MCL**, **FastTree**, **FastME**,
  and **MAFFT**. (OrthoFinder runs DIAMOND, MCL, FastTree, MAFFT, and FastME
  internally.)

**For the time-calibrated species tree:**

- **r8s 1.81** — *not available on conda*. Build it once from source with the
  helper in [`tools/r8s/`](tools/r8s/README.md). If r8s is not installed, the
  OrthoFinder job still completes normally and simply skips the dating step.

**For turnover modeling and trait association:**

- **CAFE-5** — installed separately (see the CAFE-5 project). Used by the
  turnover-modeling stage.
- **R** for the downstream R steps — **included in the `convgeno` conda env**
  (`environment.yml` adds `r-base` + **ape**, **castor**, **expm**, **yaml**), so a
  single `conda activate convgeno` covers OrthoFinder *and* the automatic
  categorical phenotype tree (CAFE `-y`) step end-to-end. Its ancestral-state
  reconstruction is vendored (`Src/Reu/r/anc_recon.R`), so **RERconverge is not
  required**. To run that step against an external/module-provided R without conda,
  install the same stack with [`tools/r-deps/`](tools/r-deps/README.md). The planned
  association stage additionally uses **RERconverge**, **ggplot2**, and **ggtree**
  (see note below).

> **RERconverge (optional here):** the vendored engine already reproduces the
> Meyer Lab reconstruction (`getAncLiks`). To run the lab's own build instead —
> or for the association stage that targets its functions — install/attach
> RERconverge (the lab fork, not the CRAN/official release) or source the lab's
> `RERConvergeFunctions.R` and set `phenotype_tree.rer_functions`.

---

## Installation

```bash
# 1. Create and activate the conda environment (installs OrthoFinder + friends)
conda env create -f environment.yml
conda activate convgeno

# 2. Editable-install the convgeno Python package (also provides the CLI)
pip install -e Src/Reu/python/
```

The environment is named **`convgeno`** and provides the `convgeno` command-line
tool. Then, once, build r8s for the dating step:

```bash
# See tools/r8s/README.md for module/toolchain details
./tools/r8s/install_r8s.sh
```

---

## Quick start

From a login node with the environment active:

```bash
# 1. Clean your proteomes: keep the longest protein per gene, one FASTA per species
python Src/Loc/scripts/filter_isoforms.py dir raw_proteomes/ cleaned_proteomes/

# 2. One-time interactive setup: detects your cluster and writes pipeline_config.yaml
convgeno init

# 3. Generate + submit the OrthoFinder workflow (multi-node is the default)
convgeno orthofinder run          # shows the scripts, then asks to submit
#   …or skip the prompt:
convgeno orthofinder run -y

# 4. Watch it
squeue -u "$USER"
```

Results land in the `orthofinder.output_dir` from your config, and — if r8s is
installed and you gave a calibration during `init` — a time-calibrated
`species_tree_ultrametric.nwk` appears alongside them.

---

## Usage

Run everything from the repository root, with the `convgeno` environment active.

### 1. Filter isoforms

OrthoFinder expects one clean FASTA per species with a
single representative protein per gene. This step reduces each proteome to its
**longest isoform per gene**. It streams gzip/bz2 inputs transparently,
auto-detects Ensembl vs NCBI headers, and has a low-memory mode for very large
proteomes.

```bash
# Batch: one FASTA per species in a directory -> a directory of filtered FASTAs
python Src/Loc/scripts/filter_isoforms.py dir <input_dir> <output_dir> \
    [--stats-json stats.json]

# Single species from one or more files (e.g. per-chromosome) -> one FASTA
python Src/Loc/scripts/filter_isoforms.py files -i in1.fa[.gz] [in2.fa …] -o out.fa
```

Shared flags:

| Flag | Meaning |
|---|---|
| `--format {auto,ensembl,ncbi}` | Header style to parse (default `auto`). |
| `--on-duplicate {error,warn,skip}` | What to do when two records share a gene ID. |
| `--stats-json PATH` | Write per-species kept/dropped counts (batch mode). |
| `--log-level`, `--log-file` | Logging verbosity and destination. |

Tip labels downstream come from the **file basename** (without extension), so
name each output file for its species — e.g. `Homo_sapiens.fa`. Recognized
extensions are `.fa`, `.fasta`, and `.faa`.

### 2. Set up your run (`convgeno init`)

`convgeno init` is a one-time interactive wizard. It probes
your cluster (partitions, cores, memory, scratch) and records the absolute conda
paths that generated SLURM scripts need to activate the environment on bare
compute nodes. It writes everything to **`pipeline_config.yaml`**.

```bash
convgeno init                 # writes ./pipeline_config.yaml
convgeno init --output my_config.yaml
```

During setup you also choose how the **species tree** is supplied and (optionally)
give one **divergence-time calibration** so the tree can be dated in real time:

- **Species tree** — let OrthoFinder build it (default), or bring your own
  (ultrametric or not). See
  [`docs/species_tree/species_tree.md`](docs/species_tree/species_tree.md) for
  the three options and the exact flags.
- **Calibration** — name two species in your data and the divergence time (in
  millions of years) between them; r8s uses this one number to scale the whole
  tree. You can skip it (you get a relative-time tree) or add several
  calibrations later in the config. The full walkthrough is in
  [`docs/r8s/divergence-time-calibration.md`](docs/r8s/divergence-time-calibration.md).

### 3. Run OrthoFinder and build the time-calibrated tree

`convgeno orthofinder` generates and submits the SLURM
job(s) that run OrthoFinder. **Multi-node is the default**; single-node is
opt-in for benchmarking or small runs.

```bash
convgeno orthofinder generate          # write the script(s), do NOT submit
convgeno orthofinder run               # write + submit (prompts for confirmation)
convgeno orthofinder run -y            # …skip the prompt
convgeno orthofinder run --single-node # single-node benchmark / fallback
```

Options:

| Flag | Meaning |
|---|---|
| `--mode {multinode,single-node}` | Execution mode (default `multinode`). |
| `--multinode` / `--single-node` | Aliases for `--mode`. Conflicting flags error out. |
| `--config PATH` | Config file to read (default: the mode-specific `pipeline_config_<mode>.yaml` if present, else `pipeline_config.yaml`). |
| `--script-dir DIR` | Where multi-node scripts are written (default `slurm_scripts/`). |
| `--script PATH` | Output path for the single-node script (default `slurm_scripts/orthofinder.sh`). |
| `-y`, `--yes` | Submit without the confirmation prompt (`run` only). |

**What you'll see when it runs.** OrthoFinder has no built-in multi-node mode,
so the multi-node path splits the run into **three chained SLURM jobs** using
OrthoFinder's own checkpoints — prepare (`-op`) → a search array
(`diamond blastp`) → resume (`-b`):

```
Job A: prepare  ──afterok──▶  Job B: search array  ──afterok──▶  Job C: resume
```

They appear together in `squeue`. The search work is spread across array tasks
with a longest-processing-time balancer (Graham, 1969), and the final resume
step uses the **same gene-tree method as a single-node run**, so the orthogroups
and trees are identical either way. If your cluster rejects the SLURM account in
your config, submission automatically retries without it. Phase-by-phase detail:
[prepare](docs/orthofinder/1-prepare-phase.md) ·
[search](docs/orthofinder/2-search-phase.md) ·
[resume](docs/orthofinder/3-resume-phase.md).

**The dated tree comes for free.** At the tail of a successful OrthoFinder job,
the pipeline runs **r8s** to turn OrthoFinder's species tree (branch lengths in
substitutions per site) into a rooted, ultrametric **time tree** — the input
CAFE-5 needs — using the calibration you set in `init`. The result,
`species_tree_ultrametric.nwk`, is written into the OrthoFinder output
directory, along with the exact r8s control file and log for inspection. This
step never risks your OrthoFinder results: it runs only after OrthoFinder has
already succeeded, and if r8s is missing it is skipped with a clear message. See
[`docs/r8s/divergence-time-calibration.md`](docs/r8s/divergence-time-calibration.md).

> Advanced users can run the dating step by hand with
> `Src/Loc/scripts/make_tree_ultrametric.py` (`-p 'sp1,sp2' -c <Myr>`, or
> repeatable `--calibration NAME:SP1,SP2:AGE`, plus `--nsites`, `--smoothing`,
> `--cross-validate`, `--r8s-path`, `--dry-run`). See the same doc.

### 4. Model gene-family turnover with CAFE-5

CAFE-5 takes the gene-family count matrix and the time-calibrated species tree
and models where each family expands or contracts along the tree. The pipeline
builds the CAFE-5 count matrix from OrthoFinder's orthogroups
(`Src/Loc/scripts/prepare_cafe_inputs.py`) and pairs it with the
`species_tree_ultrametric.nwk` from stage 3 — which is already the **binary,
rooted, ultrametric time tree** that CAFE-5 requires. For categorical
phenotypes, a helper builds the multi-λ (`-y`) rate tree CAFE-5 uses to let
gain/loss rates differ between trait states. CAFE-5 is installed separately.

### 5. Associate copy number with convergent traits

The final stage links gene-family copy-number change to convergent phenotypes
using a **RERconverge**-based branch-association approach, then produces summary
tables and figures. Reusable R utilities for parsing OrthoFinder and CAFE-5
output, building phenotype trees, and plotting live in
[`Src/Reu/r/`](Src/Reu/r/).

### Calibrating DIAMOND throughput

The multi-node search array estimates each task's walltime from a per-core
DIAMOND rate, `throughput_const`. A default
calibrated on the Lehigh Sol cluster ships in the code, but DIAMOND speed varies
by CPU. To size walltimes accurately on **your** cluster, measure the constant
with the kit in
[`tools/throughput-constant/`](tools/throughput-constant/README.md) and set
`multinode.throughput_const` in your config. This is optional — the built-in
default plus a safety margin works out of the box.

---

## Inputs and outputs

**What you provide:**

- One protein FASTA **per species** (raw proteomes; the isoform filter cleans
  them). File basenames become the species/tip labels everywhere downstream, so
  name them meaningfully (`Homo_sapiens.fa`, `Felis_catus.fa`, …).
- *(Optional)* Your own species tree in Newick format, and/or one or more
  divergence-time calibrations — supplied during `convgeno init`.

**What you get** (in `orthofinder.output_dir`):

- OrthoFinder's full results — orthogroups (gene families), gene trees, the
  species tree, and comparative-genomics statistics.
- `species_tree_ultrametric.nwk` — the time-calibrated species tree (when r8s is
  installed), plus `r8s_work/` with the control file and r8s log.
- A gene-family count matrix suitable for CAFE-5.

---

## Configuration

`convgeno init` writes `pipeline_config.yaml`. You normally never edit it by
hand, but these are the blocks you may want to tune:

```yaml
project_dir: /path/to/project
conda_env: convgeno

slurm:                       # cluster resources (auto-detected at init)
  partition: hawkcpu
  cpus_per_task: 48
  mem: 350400M
  time_limit: "72:00:00"
  account: null              # blank / null / none -> no --account line emitted
  scratch_dir: /path/scratch # optional

orthofinder:
  input_dir: cleaned_proteomes          # one filtered FASTA per species
  output_dir: Data/processed/orthofinder_<timestamp>
  sequence_search: diamond
  msa_program: mafft         # OrthoFinder -A
  tree_program: fasttree     # OrthoFinder -T

ultrametric:                 # the r8s dating step (from your init answers)
  species_a: Homo_sapiens
  species_b: Felis_catus
  divergence_my: 94.0
  smoothing: 100             # higher = closer to a strict molecular clock

multinode:                   # OPTIONAL — every knob auto-derives if omitted
  throughput_const: 6.055387e11   # from the throughput-calibration kit
```

The `multinode:` block is entirely optional: leave it out and the search-array
sizing is derived automatically from your cluster and proteome sizes. The
`ultrametric:` block also accepts a list of `calibrations:` with fixed ages or
age windows — see the
[divergence-time doc](docs/r8s/divergence-time-calibration.md).

---

## Command reference

| Command | What it does |
|---|---|
| `convgeno init [--output PATH]` | Interactive first-time setup; writes `pipeline_config.yaml`. |
| `convgeno orthofinder generate [--mode …] [--config …] [--script-dir …]` | Write the OrthoFinder SLURM script(s) **without** submitting. |
| `convgeno orthofinder run [--mode …] [-y] [--config …] [--script-dir …]` | Write **and** submit the OrthoFinder job(s). |
| `python Src/Loc/scripts/filter_isoforms.py dir <in> <out>` | Batch longest-isoform filtering (directory of proteomes). |
| `python Src/Loc/scripts/filter_isoforms.py files -i <in…> -o <out>` | Longest-isoform filtering from one or more files into one FASTA. |
| `python Src/Loc/scripts/make_tree_ultrametric.py …` | Run the r8s dating step by hand (normally automatic). |

`convgeno` (no arguments) prints help. See [Usage](#usage) for the full flag
lists.

---

## Documentation

| Topic | Where |
|---|---|
| Multi-node OrthoFinder — prepare phase | [`docs/orthofinder/1-prepare-phase.md`](docs/orthofinder/1-prepare-phase.md) |
| Multi-node OrthoFinder — search phase | [`docs/orthofinder/2-search-phase.md`](docs/orthofinder/2-search-phase.md) |
| Multi-node OrthoFinder — resume phase | [`docs/orthofinder/3-resume-phase.md`](docs/orthofinder/3-resume-phase.md) |
| Choosing / supplying the species tree | [`docs/species_tree/species_tree.md`](docs/species_tree/species_tree.md) |
| Divergence-time calibration (r8s) | [`docs/r8s/divergence-time-calibration.md`](docs/r8s/divergence-time-calibration.md) |
| Installing r8s from source | [`tools/r8s/README.md`](tools/r8s/README.md) |
| Installing the R packages | [`tools/r-deps/README.md`](tools/r-deps/README.md) |
| Calibrating the DIAMOND throughput constant | [`tools/throughput-constant/README.md`](tools/throughput-constant/README.md) |

---

## How to cite

This pipeline is orchestration around established tools and methods. If you use
it, please cite the tools for the stages you actually ran:

**Gene-family (orthogroup) inference — OrthoFinder**

- Emms, D. M., & Kelly, S. (2019). OrthoFinder: Phylogenetic orthology inference
  for comparative genomics. *Genome Biology, 20*, 238.
  https://doi.org/10.1186/s13059-019-1832-y
- Emms, D. M., & Kelly, S. (2015). OrthoFinder: Solving fundamental biases in
  whole-genome comparisons dramatically improves orthogroup inference accuracy.
  *Genome Biology, 16*, 157. https://doi.org/10.1186/s13059-015-0721-2

  *OrthoFinder internally runs DIAMOND, MCL, FastTree, MAFFT, and FastME
  (installed by the conda environment); cite these as recommended in the
  OrthoFinder documentation if you use those components.*

**Time-calibrated species tree — r8s**

- Sanderson, M. J. (2003). r8s: Inferring absolute rates of molecular evolution
  and divergence times in the absence of a molecular clock. *Bioinformatics,
  19*(2), 301–302. https://doi.org/10.1093/bioinformatics/19.2.301

**Gene-family turnover modeling — CAFE-5**

- Mendes, F. K., Vanderpool, D., Fulton, B., & Hahn, M. W. (2020). CAFE 5 models
  variation in evolutionary rates among gene families. *Bioinformatics,
  36*(22–23), 5516–5518. https://doi.org/10.1093/bioinformatics/btaa1022

**Convergent-trait association — RERconverge**

- Kowalczyk, A., Meyer, W. K., Partha, R., Mao, W., Clark, N. L., & Chikina, M.
  (2019). RERconverge: An R package for associating evolutionary rates with
  convergent traits. *Bioinformatics, 35*(22), 4815–4817.
  https://doi.org/10.1093/bioinformatics/btz468
- Partha, R., Kowalczyk, A., Clark, N. L., & Chikina, M. (2019). Robust methods
  for detecting convergent shifts in evolutionary rates. *Molecular Biology and
  Evolution, 36*(8), 1817–1830. https://doi.org/10.1093/molbev/msz107
- Redlich, R., Kowalczyk, A., Tene, M., Sestili, H. H., Foley, K., Saputra, E.,
  Clark, N. L., Chikina, M., Meyer, W. K., & Pfenning, A. R. (2024). RERconverge
  expansion: Using relative evolutionary rates to study complex categorical
  trait evolution. *Molecular Biology and Evolution, 41*(11), msae210.
  https://doi.org/10.1093/molbev/msae210

**Downstream phylogenetics & figures (R)**

- Paradis, E., & Schliep, K. (2019). ape 5.0: An environment for modern
  phylogenetics and evolutionary analyses in R. *Bioinformatics, 35*(3),
  526–528. https://doi.org/10.1093/bioinformatics/bty633
- Louca, S., & Doebeli, M. (2018). Efficient comparative phylogenetics on large
  trees (castor). *Bioinformatics, 34*(6), 1053–1055.
  https://doi.org/10.1093/bioinformatics/btx701
- Wickham, H. (2016). *ggplot2: Elegant graphics for data analysis*. Springer.
  https://doi.org/10.1007/978-3-319-24277-4
- Yu, G., Smith, D. K., Zhu, H., Guan, Y., & Lam, T. T.-Y. (2017). ggtree: An R
  package for visualization and annotation of phylogenetic trees with their
  covariates and other associated data. *Methods in Ecology and Evolution,
  8*(1), 28–36. https://doi.org/10.1111/2041-210X.12628

**Multi-node search scheduling (this pipeline's load balancer)**

- Graham, R. L. (1969). Bounds on multiprocessing timing anomalies. *SIAM
  Journal on Applied Mathematics, 17*(2), 416–429.
  https://doi.org/10.1137/0117039
