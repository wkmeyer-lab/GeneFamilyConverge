# GeneFamilyConverge

**A pipeline for characterizing gene-family copy number across species and
associating copy-number variation with convergent phenotypes.**

Supports **Aim 1a** of the Meyer Lab NSF CAREER project. The implemented core is
the `convgeno` Python package: longest-isoform proteome filtering and a
**multi-node, SLURM-orchestrated OrthoFinder workflow** that distributes the
O(n²) all-vs-all sequence search across a cluster while keeping results
scientifically equivalent to a single-node run.

---

## Scientific overview

```
proteomes ──▶ filter isoforms ──▶ OrthoFinder ──▶ gene-family counts
 (1/species)   (longest per gene)  (orthogroups)        │
                                                        ▼
                                                     CAFE-5
                                                 (size/turnover)
                                                        │
                                                        ▼
                                            R analysis + figures
                                     (copy-number ↔ convergent phenotype)
```

The scientific question is whether **shifts in gene-family copy number**
coincide with **independently evolved (convergent) phenotypes**. OrthoFinder
defines the gene families (orthogroups); a phylogenetic turnover model (CAFE-5)
characterizes expansions/contractions along a species tree; and a downstream R
step links copy-number change to convergent traits (the planned method is a
RERconverge-adapted branch association).

## Implementation status

This repository is **partially implemented** — the honest picture:

| Stage | Status | Where |
|---|---|---|
| Longest-isoform filtering | ✅ **Implemented & tested** | `convgeno.io.fasta`, `Src/Loc/scripts/filter_isoforms.py` |
| OrthoFinder — multi-node SLURM orchestration (prepare → search → resume) | ✅ **Implemented & tested** | `convgeno.slurm`, `convgeno orthofinder` |
| OrthoFinder — single-node (benchmark/fallback) | ✅ **Implemented & tested** | `convgeno orthofinder … --single-node` |
| Cluster auto-discovery + throughput calibration | ✅ **Implemented** | `convgeno.slurm.discovery`, `tools/throughput-constant/` |
| CAFE-5 turnover modeling | ⏳ **Planned (stubs)** | `convgeno.external.cafe` |
| R association analysis + figures | ⏳ **Planned (stubs)** | `Src/Reu/r/`, `Src/Loc/scripts` |

Everything below documents the **implemented** parts.

---

## Requirements

- A **SLURM** HPC cluster (for the OrthoFinder workflow) with `sinfo`,
  `scontrol`, `sbatch` (and optionally `sacctmgr`).
- **conda** (miniforge/miniconda).
- Python **3.10–3.12**.
- Bioinformatics tools (installed by the conda env): **OrthoFinder 2.5.5**,
  **DIAMOND 2.1.9**, **MCL**, **FastTree**, **FastME**, **MAFFT**.
- **r8s 1.81** for the ultrametric step — **not** on conda; build it from source
  with [`tools/r8s/install_r8s.sh`](tools/r8s/README.md).

## Installation

```bash
conda env create -f environment.yml
conda activate convgeno
pip install -e Src/Reu/python/     # editable install of the convgeno package
```

The env is named **`convgeno`** and installs the tools above plus the Python
package (with its `convgeno` CLI).

> **r8s** (needed to make the species tree ultrametric for CAFE-5) is not on
> conda. After the env is set up, build it once with
> `./tools/r8s/install_r8s.sh` — see [`tools/r8s/README.md`](tools/r8s/README.md).

> On a cluster, run **`convgeno init`** once (below). It detects your partition,
> cores, memory, QOS, and scratch, and records absolute conda paths so generated
> SLURM scripts can bootstrap the environment on bare compute nodes.

---

## Quick start (multi-node OrthoFinder)

From a directory containing your proteomes, on a login node with the env active:

```bash
# 1. One-time interactive setup: detects the cluster, writes pipeline_config.yaml
convgeno init

# 2. Generate + submit the whole chain (multi-node is the DEFAULT).
#    Submits prepare -> search array -> resume, chained by afterok.
convgeno orthofinder run          # shows the 3 scripts, then prompts
#   …or non-interactive:
convgeno orthofinder run -y

# 3. Watch it
squeue -u "$USER"
```

Results land in the configured `orthofinder.output_dir` (`Results_*/`). The
single-node benchmark path is `convgeno orthofinder run --single-node -y`.

---

## Usage

### 1. Longest-isoform filtering

Reduces each proteome to one representative (longest) sequence per gene — the
only *scientific* preprocessing step, and the input OrthoFinder expects (one
clean FASTA per species). It streams gzip/bz2 transparently (magic-byte
sniffing), auto-detects Ensembl vs NCBI headers, and offers a low-memory
two-pass mode. It is a `Src/Loc` script that imports `convgeno.io.fasta`:

```bash
# Batch: one FASTA per species in a directory -> filtered dir (+ optional stats)
python Src/Loc/scripts/filter_isoforms.py dir <input_dir> <output_dir> \
    [--stats-json stats.json]

# Single species from one or more inputs (e.g. per-chromosome) -> one output
python Src/Loc/scripts/filter_isoforms.py files -i in1.fa[.gz] [in2.fa …] -o out.fa

# Shared flags: --format {auto,ensembl,ncbi}  --on-duplicate {error,warn,skip}
#               --log-level  --log-file
```

### 2. OrthoFinder (multi-node — default)

`convgeno orthofinder` generates and submits SLURM jobs. Multi-node is the
default; single-node is opt-in.

```bash
convgeno orthofinder generate            # write the 3 scripts, do NOT submit
convgeno orthofinder run                 # write + submit the chain (prompts)
convgeno orthofinder run -y              # …skip the prompt
convgeno orthofinder run --single-node   # single-node benchmark/fallback
```

- **Mode selection:** `--mode {multinode,single-node}` (default `multinode`),
  with `--multinode` / `--single-node` as aliases. Conflicting flags error out.
- **Scripts** are written to `slurm_scripts/` (`--script-dir` to change);
  config is `pipeline_config.yaml` (`--config` to change).
- The three jobs are chained `prepare → search → resume` with `afterok`, so a
  failed link halts the chain.

See **[How the multi-node workflow works](#how-the-multi-node-workflow-works)**
and the detailed phase docs in [`docs/orthofinder/`](docs/orthofinder/).

### 3. Calibrating DIAMOND throughput

The search array sizes each task's walltime from a per-core DIAMOND rate
(`throughput_const`). A default calibrated on Lehigh Sol (`hawkcpu`) ships in the
code; recalibrate on your own cluster with the kit in
[`tools/throughput-constant/`](tools/throughput-constant/) (run the experiment,
then set `multinode.throughput_const` in your config).

---

## Configuration

`convgeno init` writes `pipeline_config.yaml`. Key blocks:

```yaml
project_dir: /path/to/project
conda_env: convgeno

slurm:                       # cluster resources (auto-derived at init)
  partition: hawkcpu
  cpus_per_task: 48          # physical cores − 4 (fat-node profile)
  mem: 350400M               # cpus × mem-per-cpu
  time_limit: "72:00:00"
  scratch_dir: /path/scratch # optional; enables scratch for search + resume temp
  is_ephemeral_scratch: false

orthofinder:
  input_dir: Data/interim/cleaned_proteomes   # one FASTA per species
  output_dir: Data/processed/orthofinder_multinode_<timestamp>
  sequence_search: diamond
  msa_program: mafft         # -M msa -A mafft
  tree_program: fasttree     # -T fasttree
  search_threads: 48         # -t
  analysis_threads: 12       # -a

runtime:                     # absolute conda paths for compute-node bootstrap
  conda_module: …
  conda_base: …
  conda_env_prefix: …

multinode:                   # OPTIONAL — all auto-derived if omitted
  # search_cpus, array_throttle, waves, threads_per_command,
  # throughput_const, time_margin, search_time_limit, search_mem
  throughput_const: 6.055387e11   # override the built-in calibrated default
```

The `multinode:` block is entirely optional: leave it out and every search-array
knob is auto-derived from cluster discovery + proteome sizes. Any field pins that
value.

---

## How the multi-node workflow works

OrthoFinder has no built-in multi-node mode, but exposes checkpoints — `-op`
(prepare) and `-b` (resume) — that isolate the one embarrassingly-parallel phase
(the search). `convgeno` splits the run into three chained SLURM jobs:

```
Job 1 PREPARE ──afterok──▶ Job 2 SEARCH (array) ──afterok──▶ Job 3 RESUME
 (1 node, orthofinder -op)   (many tasks, diamond blastp × n²)  (1 fat node, -b)
```

- **Prepare** runs `orthofinder -op`, builds the per-species DIAMOND databases,
  asserts completeness (exactly `n²` searches, `n` databases), and writes a
  **cost-balanced per-task manifest**.
- **Search** is a SLURM array of `T` LPT-balanced tasks
  (`cost(i,j)=|Sᵢ|·|Sⱼ|`), each running `C/p` `diamond blastp` concurrently, with
  idempotent skip of completed results and optional scratch staging.
- **Resume** verifies all `n²` results exist, raises the open-file limit to the
  `≈n²` requirement (OrthoFinder issue #571), then runs `orthofinder -b` with the
  **same `-M/-A/-T` gene-tree method as single-node** — guaranteeing identical
  orthogroups and trees.

Full detail, math, and runtime behavior:

- [`docs/orthofinder/1-prepare-phase.md`](docs/orthofinder/1-prepare-phase.md)
- [`docs/orthofinder/2-search-phase.md`](docs/orthofinder/2-search-phase.md)
- [`docs/orthofinder/3-resume-phase.md`](docs/orthofinder/3-resume-phase.md)

All parallelism is derived generically from **species count, per-species
proteome sizes, and live cluster discovery** — nothing is hardcoded or tuned from
a prior run.

---

## Repository layout

```
Src/Reu/     Reusable library code — NO hardcoded paths, everything by argument
  python/    the installable `convgeno` package (+ tests):
    convgeno/io/fasta.py        longest-isoform filtering
    convgeno/slurm/             cluster discovery, config, SLURM generators
      multinode_generator.py    prepare / search-array / resume generators + sizing
      script_generator.py       single-node OrthoFinder generator
      discovery.py              sinfo/scontrol/sacctmgr probing
      config.py                 SlurmConfig / PipelineConfig / MultinodeConfig
      build_search_manifest.py  prepare-time LPT manifest builder (python -m)
      verify_search_complete.py resume completeness gate (python -m)
    convgeno/cli/               argparse entry point + init / orthofinder handlers
  r/         standalone R utilities (planned)
Src/Loc/     Project-specific glue (this run's scripts/configs; allowed to hardcode)
Data/        raw → interim → processed, plus logs/ and manifests/
docs/        documentation (docs/orthofinder/ = the phase docs)
tools/       operator tooling (throughput-constant/ = DIAMOND calibration kit)
```

**Reu vs Loc** is the core convention: `Reu` is reusable and path-free; `Loc` is
this project's glue and may reference concrete paths. Python is upstream
(filtering, format prep, tool invocation); R is downstream (parsing, stats,
figures); the two communicate **only via files** (TSV/FASTA/Newick).

---

## Testing

```bash
cd Src/Reu/python
pytest                        # full suite (550+ tests; config in pyproject.toml)
pytest --cov=convgeno         # with coverage
ruff check .                  # lint (E,F,W,I,N,UP,B,SIM; line-length 88)
ruff format .
```

Tests use `tmp_path` fixtures and never write outside their temp dir. Generated
SLURM scripts are checked with `bash -n`; the LPT balancer, sizing math,
completeness/fd gates, and CLI dispatch are unit-tested.

---

## Underlying tools & citations

If you use this pipeline, please cite the tools it invokes (and the specific
versions in `environment.yml`):

- **OrthoFinder** — gene-family (orthogroup) inference. Emms, D. M., & Kelly, S.
  (2019). *Genome Biology*. https://doi.org/10.1186/s13059-019-1832-y
- **DIAMOND** — fast protein sequence search (the O(n²) all-vs-all step; this
  pipeline uses v2.1.9). Buchfink, B., Xie, C., & Huson, D. H. (2015).
  *Nature Methods*. https://doi.org/10.1038/nmeth.3176 (v2: Buchfink, Reuter, &
  Drost, 2021, *Nature Methods*).
- **MCL** (via TRIBE-MCL) — orthogroup clustering. Enright, A. J., Van Dongen,
  S., & Ouzounis, C. A. (2002). *Nucleic Acids Research*.
  https://doi.org/10.1093/nar/30.7.1575
- **FastTree 2** — gene/species tree inference. Price, M. N., Dehal, P. S., &
  Arkin, A. P. (2010). *PLoS ONE*. https://doi.org/10.1371/journal.pone.0009490
- **MAFFT** — multiple sequence alignment (OrthoFinder `-M msa`). See
  https://mafft.cbrc.jp/alignment/software/
- **FastME** — species-tree distance method used by OrthoFinder. See
  http://www.atgc-montpellier.fr/fastme/

> Citation records were retrieved and DOI-verified via the scite literature
> service; confirm volumes/pages against the publisher before publishing.

---

## Roadmap

Implemented: isoform filtering and the multi-node/single-node OrthoFinder
orchestration. Planned (currently stubs): CAFE-5 turnover modeling and the R
association analysis (RERconverge-adapted branch association) linking
copy-number change to convergent phenotypes.

## Acknowledgements

Developed for the **Meyer Lab** NSF CAREER project (Aim 1a), Lehigh University.
